#ifndef PUMIPIC_ADJACENCY_HPP
#define PUMIPIC_ADJACENCY_HPP

#include <iostream>

#include "Omega_h_for.hpp"
#include "Omega_h_adj.hpp"
#include "Omega_h_element.hpp"

#include "pumipic_utils.hpp"
#include "pumipic_constants.hpp"


namespace pumipic
{
  const int verbose = 0; //TODO move to pumipic_constants.hpp as inline/static

/*
   see description: Omega_h_simplex.hpp, Omega_h_refine_topology.hpp line 26
   face_vert:0,2,1; 0,1,3; 1,2,3; 2,0,3.
   corresp. opp. vertexes: 3,2,0,1, by simplex_opposite_template(DIM, FDIM, iface, i) ?
   side note: r3d.cpp line 528: 3,2,1; 0,2,3; 0,3,1; 0,1,2 .Vertexes opp.:0,1,2,3
              3
            / | \
          /   |   \
         0----|----2
          \   |   /
            \ | /
              1
*/
//retrieve face coords in the o order
OMEGA_H_INLINE void get_face_coords(const o::Matrix<DIM, 4> &M,
          const o::LO iface, o::Few<o::Vector<DIM>, 3> &abc) {

   //face_vert:0,2,1; 0,1,3; 1,2,3; 2,0,3
    OMEGA_H_CHECK(iface<4 && iface>=0);
    abc[0] = M[o::simplex_down_template(DIM, FDIM, iface, 0)];
    abc[1] = M[o::simplex_down_template(DIM, FDIM, iface, 1)];
    abc[2] = M[o::simplex_down_template(DIM, FDIM, iface, 2)];

    if(verbose > 1)
        std::cout << "face " << iface << ": \n"; 
}

OMEGA_H_INLINE void get_edge_coords(const o::Few<o::Vector<DIM>, 3> &abc,
          const o::LO iedge, o::Few<o::Vector<DIM>, 2> &ab) {

   //edge_vert:0,1; 1,2; 2,0
    ab[0] = abc[o::simplex_down_template(FDIM, 1, iedge, 0)];
    ab[1] = abc[o::simplex_down_template(FDIM, 1, iedge, 1)];
    if(verbose > 2)
        std::cout << "abc_index " << ab[0].data() << ", " << ab[1].data()
                  << " iedge:" << iedge << "\n";
}

OMEGA_H_INLINE void check_face(const o::Matrix<DIM, 4> &M,
    const o::Few<o::Vector<DIM>, 3>& face, const o::LO faceid ){

    o::Few<o::Vector<DIM>, 3> abc;
    get_face_coords( M, faceid, abc);

    if(verbose > 2) {
      print_array(abc[0].data(),3, "a");
      print_array(face[0].data(),3, "face1");
      print_array(abc[1].data(), 3, "b");
      print_array(face[1].data(), 3, "face2");
      print_array(abc[2].data(), 3, "c");
      print_array(face[2].data(), 3, "face3");
    }
    OMEGA_H_CHECK(true == compare_array(abc[0].data(), face[0].data(), DIM)); //a
    OMEGA_H_CHECK(true == compare_array(abc[1].data(), face[1].data(), DIM)); //b
    OMEGA_H_CHECK(true == compare_array(abc[2].data(), face[2].data(), DIM)); //c
}

// BC coords are not in order of its corresp. opp. vertexes. Bccoord of tet(iface, xpoint)
//TODO Warning: Check opposite_template use in this before using
OMEGA_H_INLINE bool find_barycentric_tet( const o::Matrix<DIM, 4> &Mat,
     const o::Vector<DIM> &pos, o::Write<o::Real> &bcc) {

  for(o::LO i=0; i<3; ++i) bcc[i] = -1;

  o::Real vals[4];
  o::Few<o::Vector<DIM>, 3> abc;
  for(o::LO iface=0; iface<4; ++iface) {// TODO last not needed
    get_face_coords(Mat, iface, abc);
    auto vab = abc[1] - abc[0]; //b - a;
    auto vac = abc[2] - abc[0]; //c - a;
    auto vap = pos - abc[0]; // p - a;
    vals[iface] = osh_dot(vap, o::cross(vac, vab)); //ac, ab NOTE

    if(verbose > 2) {
      std::cout << "vol: " << vals[iface] << " for points_of_this_TET:\n" ;
      print_array(abc[0].data(),3);
      print_array(abc[1].data(),3);
      print_array(abc[2].data(),3);
      print_array(pos.data(),3, "point");
      std::cout << "\n";
    }
  }
  //volume using bottom face=0
  get_face_coords(Mat, 0, abc);
  auto vtx3 = o::simplex_opposite_template(DIM, FDIM, 0);
  OMEGA_H_CHECK(3 == vtx3);
  // abc in order, for bottom face: M[0], M[2](=abc[1]), M[1](=abc[2])
  o::Vector<DIM> cross_ac_ab = o::cross(abc[2]-abc[0], abc[1]-abc[0]); //NOTE
  o::Real vol6 = osh_dot(Mat[vtx3]-Mat[0], cross_ac_ab);
  o::Real inv_vol = 0.0;
  if(vol6 > EPSILON) // TODO tolerance
    inv_vol = 1.0/vol6;
  else {
    if(verbose > 0)  
      std::cout << "Error: Volume " << vol6 << " too small \n";
    return 0;
  }
  bcc[0] = inv_vol * vals[0]; //for face0, cooresp. to its opp. vtx.
  bcc[1] = inv_vol * vals[1];
  bcc[2] = inv_vol * vals[2];
  bcc[3] = inv_vol * vals[3]; // 1-others

  return 1; //success
}


// BC coords are not in order of its corresp. vertexes. Bccoord of triangle (iedge, xpoint)
// corresp. to vertex obtained from simplex_opposite_template(FDIM, 1, iedge) ?
OMEGA_H_INLINE bool find_barycentric_tri_simple(const o::Few<o::Vector<DIM>, 3> &abc,
     const o::Vector<3> &xpoint, o::Write<o::Real> &bc) {

  o::Vector<DIM> a = abc[0];
  o::Vector<DIM> b = abc[1];
  o::Vector<DIM> c = abc[2];
  o::Vector<DIM> cross = 1/2.0 * o::cross(b-a, c-a); //NOTE order
  o::Vector<DIM> norm = o::normalize(cross);
  o::Real area = osh_dot(norm, cross);

  if(std::abs(area) < 1e-6)  //TODO
    return 0;
  o::Real fac = 1/(area*2.0);
  bc[0] = fac * osh_dot(norm, o::cross(b-a, xpoint-a));
  bc[1] = fac * osh_dot(norm, o::cross(c-b, xpoint-b));
  bc[2] = fac * osh_dot(norm, o::cross(xpoint-a, c-a));

  return 1;
}

OMEGA_H_INLINE bool line_triangle_intx_simple(const o::Few<o::Vector<DIM>, 3> &abc,
    const o::Vector<DIM> &origin, const o::Vector<DIM> &dest,
    o::Vector<DIM> &xpoint, o::LO &edge, bool reverse=false ) {

  edge = -1;
  xpoint = {0, 0, 0};

  if(verbose > 1) {
    print_osh_vector(origin, "origin", false);
    print_osh_vector(dest, "dest");
  }
    
  //Boundary exclusion. Don't set it globally and change randomnly.
  const o::Real bound_intol = 0;//SURFACE_EXCLUDE; //TODO optimum value ?

  bool found = false;
  const o::Vector<DIM> line = dest - origin;
  const o::Vector<DIM> edge0 = abc[1] - abc[0];
  const o::Vector<DIM> edge1 = abc[2] - abc[0];
  o::Vector<DIM> normv = o::cross(edge0, edge1);

  if(reverse) {
    normv = -1*normv;
    if(verbose > 0)
      std::cout << "Surface normal reversed \n";

  }
  const o::Vector<DIM> snorm_unit = o::normalize(normv);
  const o::Real dist2plane = osh_dot(abc[0] - origin, snorm_unit);
  const o::Real proj_lined =  osh_dot(line, snorm_unit);
  const o::Vector<DIM> surf2dest = dest - abc[0];

  if(std::abs(proj_lined) >0) {
    const o::Real par_t = dist2plane/proj_lined;
    if(verbose > 2)
      std::cout << " abs(proj_lined)>0;  par_t= " << par_t << " dist2plane= "
             <<  dist2plane << "; proj_lined= " << proj_lined << ";\n";
    if (par_t > bound_intol && par_t <= 1.0) {//TODO test tol value
      xpoint = origin + par_t * line;
      o::Write<o::Real> bcc{3,0};
      bool res = find_barycentric_tri_simple(abc, xpoint, bcc);
      if(verbose > 2)
        print_array(bcc.data(), 3, "BCC");
      if(res) {
        if(bcc[0] < 0 || bcc[2] < 0 || bcc[0]+bcc[2] > 1.0) { //TODO all zeros ?
          edge = min_index(bcc.data(), 3, EPSILON); //TODO test tolerance
        }
        else {
          const o::Real proj = osh_dot(snorm_unit, surf2dest);
          if(proj >0) found = true;
          else if (proj<0) {
            if(verbose > 1)
              std::cout << "Particle Entering domain\n";
          }
          else if(almost_equal(proj,0.0)) {//TODO use tol
            if(verbose > 1)
              std::cout << "Particle path on surface\n";
          }
        }
      }
      if(verbose > 2)
        print_array(bcc.data(), 3, "BCCtri");
    }
    else if(par_t >1.0) {
      if(verbose > 0)
        std::cout << "Error** Line origin and destination are on the same side of face \n";
    }
    else if(par_t < bound_intol) {// dist2plane ~0. Line contained in plane, no intersection?
      if(verbose > 1)
        std::cout << "No/Self-intersection of ptcl origin with plane at origin. t= " << par_t << " "
                << dist2plane << " " << proj_lined << "\n";
    }
  }
  else if(verbose > 1) {
      std::cout << "Line and plane are parallel \n";
  }
  return found;
}



OMEGA_H_INLINE bool search_mesh(const Omega_h::Write<Omega_h::LO> pids, Omega_h::LO nelems, const Omega_h::Write<Omega_h::Real> &x0,
 const Omega_h::Write<Omega_h::Real> &y0, const Omega_h::Write<Omega_h::Real> &z0, 
 const Omega_h::Write<Omega_h::Real> &x, const Omega_h::Write<Omega_h::Real> &y, 
 const Omega_h::Write<Omega_h::Real> &z, const Omega_h::Adj &dual, const Omega_h::Adj &down_r2f,
 const Omega_h::Read<Omega_h::I8> &side_is_exposed, const Omega_h::LOs &mesh2verts, 
 const Omega_h::Reals &coords, const Omega_h::LOs &face_verts, Omega_h::Write<Omega_h::LO> &part_flags,
 Omega_h::Write<Omega_h::LO> &elem_ids, Omega_h::Write<Omega_h::LO> &coll_adj_face_ids, 
 Omega_h::Write<Omega_h::Real> &bccs, Omega_h::Write<Omega_h::Real> &xpoints, Omega_h::LO &loops, 
 Omega_h::LO limit=0)
{
  const auto down_r2fs = &down_r2f.ab2b;
  const auto dual_faces = &dual.ab2b;
  const auto dual_elems = &dual.a2ab;

  const int debug = 0;

  const int totNumPtcls = elem_ids.size();
  Omega_h::Write<Omega_h::LO> elem_ids_next(totNumPtcls,-1);

  //particle search: adjacency + boundary crossing
  auto search_ptcl = OMEGA_H_LAMBDA( Omega_h::LO ielem)
  {
    // NOTE ielem is taken as sequential from 0 ... is it elementID ? TODO verify it
    const auto tetv2v = Omega_h::gather_verts<4>(mesh2verts, ielem);
    const auto M = Omega_h::gather_vectors<4, 3>(coords, tetv2v);

    // parallel_for loop for groups of remaining particles in this element
    //......

    // Each group of particles inside the parallel_for.
    // TODO Change ntpcl, ip start and limit. Update global(?) indices inside.
    for(Omega_h::LO ip = 0; ip < totNumPtcls; ++ip) //HACK - each element checks all particles
    {
      //skip if the particle is not in this element or has been found
      if(elem_ids[ip] != ielem || part_flags[ip] <= 0) continue;

      if(debug)
        std::cerr << "Elem " << ielem << " ptcl:" << ip << "\n";
        
      const Omega_h::Vector<3> orig{x0[ip], y0[ip], z0[ip]};
      const Omega_h::Vector<3> dest{x[ip], y[ip], z[ip]};
      
      Omega_h::Write<Omega_h::Real> bcc(4, -1.0);

      //TESTING. Check particle origin containment in current element
      find_barycentric_tet(M, orig, bcc);
      if(debug>3 && !(all_positive(bcc.data(), 4)))
          std::cerr << "ORIGIN ********NOT in elemet_id " << ielem << " \n";
      find_barycentric_tet(M, dest, bcc);

      //check if the destination is in this element
      if(all_positive(bcc.data(), 4, 0)) //SURFACE_EXCLUDE)) TODO
      {
        // TODO interpolate Fields to ptcl position, and store them, for push
        // interpolateFields(bcc, ptcls);
        elem_ids_next[ip] = elem_ids[ip];
        part_flags.data()[ip] = -1;
        if(debug) 
        {
            std::cerr << "********found in " << ielem << " \n";
            print_matrix(M);
        }
        continue;
      }
       //get element ID
      //TODO get map from omega methods. //2,3 nodes of faces. 0,2,1; 0,1,3; 1,2,3; 2,0,3
      Omega_h::LOs fmap{2,1,1,3,2,3,0,3}; 
      auto dface_ind = (*dual_elems)[ielem];
      const auto beg_face = ielem *4;
      const auto end_face = beg_face +4;
      Omega_h::LO f_index = 0;
      bool inverse;

      for(auto iface = beg_face; iface < end_face; ++iface) //not 0..3
      {
        const auto face_id = (*down_r2fs)[iface];
        if(debug >1)  
          std::cout << " \nFace: " << face_id << " dface_ind " <<  dface_ind << "\n";

        Omega_h::Vector<3> xpoint{0,0,0};
        auto fv2v = Omega_h::gather_verts<3>(face_verts, face_id); //Few<LO, 3>

        const auto face = Omega_h::gather_vectors<3, 3>(coords, fv2v);
        Omega_h::LO matInd1 = fmap[f_index*2], matInd2 = fmap[f_index*2+1];

        if(debug >3) {
          std::cout << "Face_local_index "<< fv2v.data()[0] << " " << fv2v.data()[1] << " " << fv2v.data()[2] << "\n";
          std::cout << "Mat index "<< tetv2v[matInd1] << " " << tetv2v[matInd2] << " " <<  matInd1 << " " << matInd2 << " \n";
          std::cout << "Mat dat ind " <<  tetv2v.data()[0] << " " << tetv2v.data()[1] << " "
                   << tetv2v.data()[2] << " " << tetv2v.data()[3] << "\n";
        }


        if(fv2v.data()[1] == tetv2v[matInd1] && fv2v.data()[2] == tetv2v[matInd2])
          inverse = false;
        else // if(fv2v.data()[1] == tetv2v[matInd2] && fv2v.data()[2] == tetv2v[matInd1])
        {
          inverse = true;
        }

        //TODO not useful
        auto fcoords = Omega_h::gather_vectors<3, 3>(coords, fv2v);
        auto base = Omega_h::simplex_basis<3, 2>(fcoords); //edgres = Matrix<2,3>
        auto snormal = Omega_h::normalize(Omega_h::cross(base[0], base[1]));

        Omega_h::LO dummy = -1;
        bool detected = line_triangle_intx_simple(face, orig, dest, xpoint, dummy, inverse);
        if(debug && detected)
            std::cout << " Detected: For el=" << ielem << "\n";

        if(detected && side_is_exposed[face_id])
        {
           part_flags.data()[ip] = -1;
           for(Omega_h::LO i=0; i<3; ++i)xpoints[ip*3+i] = xpoint.data()[i];
           //store current face_id and element_ids

           if(debug)
             print_osh_vector(xpoint, "COLLISION POINT");

           elem_ids_next[ip] = -1;
           break;
         }
         else if(detected && !side_is_exposed[face_id])
         {
          //OMEGA_H_CHECK(side2side_elems[side + 1] - side2side_elems[side] == 2);
           auto adj_elem  = (*dual_faces)[dface_ind];
           if(debug)
             std::cout << "Deletected For el=" << ielem << " ;face_id=" << (*down_r2fs)[iface]
                     << " ;ADJ elem= " << adj_elem << "\n";

           elem_ids_next[ip] = adj_elem;
           break;
         }

         if(!side_is_exposed[face_id])//TODO for DEBUG
         {
           if(debug)
             std::cout << "adj_element_across_this_face " << (*dual_faces)[dface_ind] << "\n";
           const Omega_h::LO min_ind = min_index(bcc.data(), 4);
           if(f_index == min_ind)
           {
             if(debug)
               std::cout << "Min_bcc el|face_id=" << ielem << "," << (*down_r2fs)[iface]
                     << " :unused adj_elem= " << (*dual_faces)[dface_ind] << "\n";
            if(!detected)
            {
              elem_ids_next[ip] = (*dual_faces)[dface_ind];
              if(debug)
                std::cout << "...  adj_elem=" << elem_ids[ip]  <<  "\n";
            }
           }

         }

         if( !side_is_exposed[face_id])
           ++dface_ind;

         ++f_index;
      } //iface 
 
    }//ip
  };

  bool found = false;
  loops = 0;
  while(!found)
  {
    if(debug) fprintf(stderr, "------------ %d ------------\n", loops);
    //TODO check if particle is on boundary and remove from list if so.

    // Searching all elements. TODO exclude those done ?
    Omega_h::parallel_for(nelems,  search_ptcl, "search_ptcl");
    found = true;
    auto cp_elm_ids = OMEGA_H_LAMBDA( Omega_h::LO i) {
      elem_ids[i] = elem_ids_next[i];
    };
    Omega_h::parallel_for(elem_ids.size(), cp_elm_ids, "copy_elem_ids");

    //FIXME using part_flags in host is ERROR. TODO memcopy part_flags
    deviceToHostFp
    // TODO use parallel reduce, else this could be a sequential bottle-neck
    for(int i=0; i<totNumPtcls; ++i){ if(part_flags[i] > 0) {found = false; break;} }
    //Copy particle data from previous to next (adjacent) element
    ++loops;

    if(limit && loops>limit) break;
  }

  std::cerr << "search iterations " << loops << "\n";

  return found;
} //search_mes



enum class TriRegion : o::LO {
  VERTA=1, 
  VERTB,
  VERTC,
  EDGEAB,
  EDGEAC,
  EDGEBC,
  TRIFACE
};


//Ref: Real-time Collision Detection by Christer Ericson, 2005.
//ptp = ref point; ptq = nearest point on triangle; abc = triangle
OMEGA_H_INLINE o::LO findNearestPointOnTriangleOrd( const o::Few< oVector, 3> &abc, 
  const oVector &ptp, oVector &ptq) {
  oVector pta = abc[0];
  oVector ptb = abc[1];
  oVector ptc = abc[2];

  oVector vab = ptb - pta;
  oVector vac = ptc - pta;
  oVector vap = ptp - pta;
  oVector vbp = ptp - ptb;
  oVector vcp = ptp - ptc;
  o::Real d1 = osh_dot(vab, vap);
  o::Real d2 = osh_dot(vac, vap);
  o::Real d3 = osh_dot(vab, vbp);
  o::Real d4 = osh_dot(vac, vbp);
  o::Real d5 = osh_dot(vab, vcp);
  o::Real d6 = osh_dot(vac, vcp);
  o::Real vc = d1*d4 - d3*d2;
  o::Real vb = d5*d2 - d1*d6;
  o::Real va = d3*d6 - d5*d4;

  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  o::Real inv = 1 / (va + vb + vc);
  o::Real v = vb * inv;
  o::Real w = vc * inv;
  o::Real u = va * inv;
  if(v>=0 && w>=0 && v<=1 && w<=1 && u >=0 && u<=1){
    qpt =  pta + v * vab+ w * vac; 
    return TriRegion::TRIFACE;
  }
// Check if P in edge region of AB, if so return projection of P onto AB
  if (vc <= 0 && d1 >= 0 && d3 <= 0) {
    o::Real v = d1 / (d1 - d3);
    // barycentric coordinates (1-v,v,0)
    qpt = pta + v * vab; 
    return TriRegion::EDGEAB;
  }
  // Check if P in edge region of AC, if so return projection of P onto AC
  if (vb <= 0 && d2 >= 0 && d6 <= 0) {
    o::Real w = d2 / (d2 - d6);
    // barycentric coordinates (1-w,0,w)
    qpt = pta + w * vac; 
    return TriRegion::EDGEAC;
  }
  // Check if P in edge region of BC, if so return projection of P onto BC
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
    o::Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    // barycentric coordinates (0,1-w,w)
    qpt =  ptb + w * (ptc - ptb); 
    return TriRegion::EDGEBC;
  }
  // Check if P in vertex region outside A
  if (d1 <= 0 && d2 <= 0) {
    // barycentric coordinates (1,0,0)
    qpt = pta;
    return TriRegion::VERTA; 
  }
  // Check if P in vertex region outside A
  if (d3 >= 0 && d4 <= d3){ 
    // barycentric coordinates (0,1,0)
    qpt = ptb;
    return TriRegion::VERTB; 
  }
  // Check if P in vertex region outside C
  if (d6 >= 0 && d5 <= d6) { 
    // barycentric coordinates (0,0,1)
    qpt = ptc; 
    return TriRegion::VERTC;
  }
}

//Ref: Real-time Collision Detection by Christer Ericson, 2005.
//ptp = ref point; ptq = nearest point on triangle; abc = triangle
OMEGA_H_INLINE o::LO findNearestPointOnTriangle( const o::Few< oVector, 3> &abc, 
  const oVector &ptp, oVector &ptq) {
  // Check if P in vertex region outside A
  oVector pta = abc[0];
  oVector ptb = abc[1];
  oVector ptc = abc[2];

  oVector vab = ptb - pta;
  oVector vac = ptc - pta;
  oVector vap = ptp - pta;
  o::Real d1 = osh_dot(vab, vap);
  o::Real d2 = osh_dot(vac, vap);
  if (d1 <= 0 && d2 <= 0) {
    // barycentric coordinates (1,0,0)
    qpt = pta;
    return TriRegion::VERTA; 
  }

  // Check if P in vertex region outside B
  oVector vbp = ptp - ptb;
  o::Real d3 = osh_dot(vab, vbp);
  o::Real d4 = osh_dot(vac, vbp);
  if (d3 >= 0 && d4 <= d3){ 
    // barycentric coordinates (0,1,0)
    qpt = ptb;
    return TriRegion::VERTB; 
  }

  // Check if P in edge region of AB, if so return projection of P onto AB
  o::Real vc = d1*d4 - d3*d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) {
    o::Real v = d1 / (d1 - d3);
    // barycentric coordinates (1-v,v,0)
    qpt = pta + v * vab; 
    return TriRegion::EDGEAB;
  }

  // Check if P in vertex region outside C
  oVector vcp = ptp - ptc;
  o::Real d5 = osh_dot(vab, vcp);
  o::Real d6 = osh_dot(vac, vcp);
  if (d6 >= 0 && d5 <= d6) { 
    // barycentric coordinates (0,0,1)
    qpt = ptc; 
    return TriRegion::VERTC;
  }

  // Check if P in edge region of AC, if so return projection of P onto AC
  o::Real vb = d5*d2 - d1*d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) {
    o::Real w = d2 / (d2 - d6);
    // barycentric coordinates (1-w,0,w)
    qpt = pta + w * vac; 
    return TriRegion::EDGEAC;
  }

  // Check if P in edge region of BC, if so return projection of P onto BC
  o::Real va = d3*d6 - d5*d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
    o::Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    // barycentric coordinates (0,1-w,w)
    qpt =  ptb + w * (ptc - ptb); 
    return TriRegion::EDGEBC;
  }

  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  o::Real inv = 1 / (va + vb + vc);
  o::Real v = vb * inv;
  o::Real w = vc * inv;
  // u*a + v*b + w*c, u = va * inv = 1 - v - w
  qpt =  pta + v * vab+ w * vac; 
  return TriRegion::TRIFACE;
}

//gitrmMesh.hpp

/*
  Reals coords;  LOs mesh2verts;  LOs dual_faces;
  LOs dual_elems;  LOs side_is_exposed;
  LOs face_verts;  LOs down_r2f;
*/

#define MESHDATA(msh) \
  const auto coords = msh.coords(); \
  const auto mesh2verts = msh.ask_elem_verts(); \
  const auto dual_faces = msh.ask_dual().ab2b; \
  const auto dual_elems= msh.ask_dual().a2ab; \
  const auto face_verts = msh.ask_verts_of(2); \
  const auto down_r2f = msh.ask_down(3, 2).ab2b; \
  const auto side_is_exposed = mark_exposed_sides(&msh);


class gitrmMesh{
public:
  gitrmMesh(o::Mesh m): msh(m){};
  ~gitrmMesh();
  gitrMesh(const gitrMesh&) = delete;

  void initNearBdryDistData();
  void fillDistToBdryFaces();
  void convert2ReadOnly();

  o::Mesh &msh;

  // Pre-process
  o::Write<o::Real> *bdryFacesDataW;
  o::Write<o::LO> *bdryFacesDataElemsW; //TODO tag

  // Reals, LOs are const & have const cast of Kokkos View
  o::Reals *bdryFacesData;
  o::LOs *bdryFacesDataElems;
}

//gitrmMesh.cpp


void gitrMesh::initNearBdryDistData(){
    OMEGA_H_CHECK(NUM_NEAR_BDRY_ELEMS > 0 && BDRYFACE_SIZE > 0);
    bdryFacesDataW = new o::Write<o::Real>( 
      NUM_NEAR_BDRY_ELEMS*BDRYFACE_SIZE, 0, "nearBdryFacesDataW");
    OMEGA_H_CHECK(bdryFacesDataW); // NOTE 

    //TOD add tag
    msh.add_tag(o::REGION, "bdryFacesFlag", 1, o::LOs(msh.nelems(), -1));
}

gitrmMesh::~gitrmMesh(){
  if(bdryFacesDataPtr)
    delete bdryFacesDataPtr;
  if(bdryFacesData)
    delete bdryFacesData; 
}
void gitrmMesh::convert2ReadOnly() {
  OMEGA_H_CHECK(bdryFacesDataW); 
  auto bdryFacesDataPtrW = o::Write<o::LO>(NUM_NEAR_BDRY_ELEMS+1, -1, 
                      "bdryFacesDataPtrW");

  o::LOs bdryFacesDataElemsTag = mesh.get_array<o::LO>(REGION, "bdryFacesFlag");

  convert = OMEGA_H_LAMBDA(o::LO elem){
    //check tag
    auto tagE = bdryFacesDataElemsTag[elem];
    if(tagE ==-1)
      return;

    ......
    
    auto nextTagE = bdryFacesDataElemsTag[elem+1];

    for(o::LO e = tagE){
      bdryFacesData = 
    }
   };
  o::parallel_for(msh.nelems(), convert, "convert_dist_data");

  bdryFacesDataPtr = new o::Reals(bdryFacesDataPtrW);
  if(bdryFacesDataW)
    delete bdryFacesDataW;
}


OMEGA_H_INLINE auto *getBdryFacesOfElem( o:Reals &bdryFacesData,
   o:LOs &bdryFacesDataPtr, o::LO elem, o::LO &nFaces) {
  OMEGA_H_CHECK(bdryFacesData);
  OMEGA_H_CHECK(bdryFacesDataPtr);
  OMEGA_H_CHECK(BDRYFACE_SIZE > 0);
  auto start = bdryFacesDataPtr[elem];
  auto end = bdryFacesDataPtr[elem+1];
  nFaces = 1.0/BDRYFACE_SIZE * (end - start);
  auto *endP = &(bdryFacesData.data()[end]);   //FIXME ??
  OMEGA_H_CHECK(endP);
  return &(bdryFacesData.data()[start]);
}

//Order of face vertices not taken care of
OMEGA_H_INLINE void getTetFaceVectors(o::Real *fd, 0::LO fi, 
    o::Few< oVector, 3> &abc){
  o::LO start = fi*BDRYFACE_SIZE;
  o::LO end = (fi+1)*BDRYFACE_SIZE;
  OMEGA_H_CHECK(fd[start]);
  OMEGA_H_CHECK(fd[end-1]);

  for(o::LO i=0; i<3; ++i){ //Tet vertexes
    for(o::LO j=0; j<3; ++j){ //coords
    abc[i][j] = fd[start++];
  }
}


//gitrmParticles.hpp. For collecting pre-processing functions 
class gitrmPcls {

  gitrmPcls(SellCSigma<Particle>* scs): scs(scs){}
  ~gitrPcls(){delete scs;}
  SellCSigma<Particle>* scs;

  thread //???

}


OMEGA_H_INLINE o::LO gitrm_findDistanceToBdry(const gitrmMesh *gm, 
  gitrmPcls *pcl) {

  MESHDATA(gm->msh);

  const auto *bdryFacesData = gm->bdryFacesData;
  const auto *bdryFacesPtr = gm->bdryFacesPtr;
  
  auto *scs = pcl->scs;
  auto &thread = pcl->thread;

  // elem will be declared; pass scs,thread. 
  PS_PARALLEL_FOR_ELEMENTS(scs, thread, elem, {
    o::LO nFaces = 0;
    const auto *bdryFacesOfElem = getBdryFacesOfElem(bdryFacesData,
                          bdryFacesPtr, elem, nFaces);

    // pid will be declared
    PS_PARALLEL_FOR_PARTICLES(scs, thread, pid, {
      const auto ref = scs.get(pid); // FIXME 
      oVector point{0, 0, 0};
      o::Few<o::Real, nFaces> dists({0});
      o::Few< oVector, 3> face;

      for( fi = 0; fi < nFaces; ++fi ){
        getTetFaceVectors(bdryFacesOfElem, fi, face);
        findNearestPointOnTriangle(face, ref, point); 
      
        dists[face] = osh_dot(point - ref, point - ref);
      }
      o::Real dist = std::min(dists); //FIXME 
      scs.distToBdry[pid] = dist;
    });
  });
}



#define DEPTH_DIST2_BDRY 0.001 // 1mm, good for device too

OMEGA_H_INLINE void findCentroidOfTet(const Omega_h::Matrix<DIM, 4> &tet, 
  oVector &cent){
    for(o::LO i=0; i<3; ++i){
      for(o::LO j=0; j<4; ++j){
          cent[i] = tet[j][i];
      }
      cent[i] *= 1/4.0;
    }
}




OMEGA_H_INLINE bool if_exists(o::LO &ids, eid, size, entry){
  for(o::LO i=0; i<size; ++i){
    if(ids[eid*NUM_BDRY_FACES + i] == entry){
      return true;
    }
  }
  return false;
}




const o::LO NUM_BDRY_FACES = 100;
const o::LO SIZE_PER_FACE = 9; //3 vtx
// mesh data, move to gitrmMesh
mesh.add_tag(REGION, "near_bdry", 1);
auto nearBdryFlags = mesh.get_array<o::LO>(REGION, "near_bdry");

//TODO add mesh tag and get these data 
o::Write<o::LO> bdryFacesUpdate(msh.nelems(), 0);
o::Write<o::LO> numBdryFaces(msh.nelems(), 0);
o::Write<o::LO> bdryFaceIds(NUM_BDRY_FACES*msh.nelems(), 0);
//Keep separate copy of bdry face coords per elem, to avoid  accessing another
// element's face. Can save storage if common bdry face data is used by face id?
o::Write<o::Real> bdryFaceData(SIZE_PER_FACE*NUM_BDRY_FACES*msh.nelems(), 0);

// Kokkos::View<double[msh.nelems()][NUM_BDRY_FACES][SIZE_PER_FACE]> bdryFaceData("bdryFaceData");


//Not checking if id already in.
OMEGA_H_INLINE bool addFaceToBdryData(data, ids, fnums, size, dof, fi, fid, elem, face){
  OMEGA_H_CHECK(fi < fnums);
  //memcpy ?
  for(o::LO i=0; i<size; ++i){
    for(o::LO j=0; j<3; j++){
      data[elem*fnums*size + fi*size + i*dof + j] = face[i][j];
    }
  }
  ids[elem*fnums + fi] = fid;  
}

// Total exposed faces has to be passed in as nbdry; no separate checking 
OMEGA_H_INLINE bool updateAdjElemFlags(dual_elems, dual_faces, elem, 
  bdryFlags, nbdry=0){

  auto dface_ind = dual_elems[elem];
  for(o::LO i=0; i<4-nbdry; ++i){
    auto adj_elem  = dual_faces[dface_ind];
    o::LO val = 1;
    Kokkos::atomic_exchange( &bdryFlags[adj_elem], val);
    ++dface_ind;
  }
}

//Before passing over to adj. elements. NOTE: Not added up.
void gitrmMesh::copyBdryFacesToSameElem(gitrmMesh *gm) {
  MESHDATA(gm->msh);

  auto fillBdry = OMEGA_H_LAMBDA(o::LO elem){

    // Checking bdry faces
    o::LO nbdry = 0; 
    o::LO bdryFids[4];
    const auto beg_face = elem *4;
    const auto end_face = beg_face +4;
    for(o::LO fi = beg_face; fi < end_face; ++fi, ++ii){//not 0..3
      const auto fid = down_r2fs[fi];
      if(side_is_exposed[fid]) {
        bdryFids[nbdry++] = fid;
      }
    }
    //If no bdry faces to add
    if(!nbdry){
      return;
    }

    //Only boundary faces
    for(o::LO fn = 0; fn < nbdry; ++fn){
      auto fid = bdryFids[fn];
      auto fv2v = o::gather_verts<3>(face_verts, fid); //Few<LO, 3>
      const auto face = o::gather_vectors<3, 3>(coords, fv2v);
      //From self bdry
      addFaceToBdryData(bdryFaceData, bdryFaceIds, NUM_BDRY_FACES, SIZE_PER_FACE, 
        fn, fid, elem, face);
    }
    numBdryFaces[elem] = nbdry;
    //Set neigbhor's flags
    updateAdjElemFlags(dual_elems, dual_faces, elem, bdryFacesUpdate, nbdry);
  };
  o::parallel_for(msh.nelems(), fillBdry, "CopyBdryFacesToSelf");
}



OMEGA_H_INLINE bool checkIfFaceWithinDistToTet(const o::Matrix<DIM, 4> &tet, 
  const o::Write<o::Real> &data, const o::LO start, const o::Real depth){
  
  o::Few<o::Vector<3>, 3> face;
  for(o::LO i=0; i<9; ++i){
    //face[i] = {data[start+i*3], data[start+i*3+1], data[start+i*3+2]};
    for(o::LO j=0; j<9; ++j){
      face[i][j] = data[sstart+i]; 
    }
  }

  for(o::LO i=0; i<3; ++i){
    for(o::LO j=0; j<4; ++j){
      oVector dv = face[i] - tet[j];
      o::Real d2 = osh_dot(dv, dv);
      if(d2 <= depth*depth){
        return true;
      }
    }
  }
  return false;
}


// Space for a fixed # of Bdry faces is assigned per element, rather than 
// using a common data to be accessed using face id. 
// First stage, only boundary faces are added to the same elements.
// When new faces are added, the owner element updates flags of adj.elements.
// So, the above inital stage sets up flags of adj.eleemnts of bdry face elements. 
// Second stage is updating and passing these faces to second and further adj. levels.
// Each element checks if its flag is set by any adj.element. If so, check all 
// adj. elements for new faces and copies off them.
// Flags are reset before checking adj.elements such that any other thread can still
// set it during copying, in which case the next iteration will be run, even if the 
// data is already copied off in the previous step due to flag/data mismatch.
typedef Kokkos::DefaultExecutionSpace exec_space; //TODO remove if duplicate
//pre-process, for dist. to bdry calculation
void gitrmMesh::addBdryFacesWithinDepth(gitrmMesh *gm, 
    o::Real depth=DIST_BDRY_DEPTH) {
  o::Mesh &msh = gm->msh;
  MESHDATA(msh);

  copyBdryFacesToSameElem(gm);

  auto fill = OMEGA_H_LAMBDA(o::LO elem){
    // This has to be init by copyBdryFacesToSameElem.
    o:LO update = bdryFacesUpdate[elem];
    //None of neigbhors have update
    if(! update){
      return;
    }
    o::LO val = 0;
    Kokkos::atomic_exchange(&bdryFacesUpdate[elem], val);

    const auto tetv2v = Omega_h::gather_verts<4>(mesh2verts, elem);
    const auto tet = Omega_h::gather_vectors<4, 3>(coords, tetv2v);

    o::LO nbf = 0;
    const auto beg_face = elem *4;
    const auto end_face = beg_face +4;
    auto dface_ind = dual_elems[elem]; //first interior
    for(auto fi = beg_face; fi < end_face; ++fi){//not 0..3
      const auto fid = down_r2fs[fi];
      if(side_is_exposed[fid]) {
        ++nbf;
        continue;
      }

      //Check all adjacent elements
      auto adj_elem  = dual_faces[dface_ind];

      auto fv2v = o::gather_verts<3>(face_verts, fid); //Few<LO, 3>
      const auto face = o::gather_vectors<3, 3>(coords, fv2v);

      bool found = false;
      o::LO sizeThis = numBdryFaces[elem];
      o::LO sizeAdj = numBdryFaces[adj_elem]; 
      o::LO tot = sizeThis;
      for(o::LO ai=0; ai < sizeAdj; ++ai){
        auto adj_fid = bdryFaceIds[adj_elem+ai]; 
        //Better way ?
        for(auto fi = 0; fi < sizeThis; ++fi){
          if(adj_fid == bdryFaceIds[elem+fi]){ 
            found = true;
            break;
          }
        }
        if(found) { 
          continue;
        }
        // In current element
        o::LO start = elem*NUM_BDRY_FACES*SIZE_PER_FACE + ai*SIZE_PER_FACE;
        bool add = checkFaceIfWithinDistToTet(tet, bdryFaceData, start, depth);
        if(!add) { 
          continue;
        }       

        addFaceToBdryData(bdryFaceData, bdryFaceIds, NUM_BDRY_FACES, SIZE_PER_FACE, 
          tot, fid, elem, face);
        ++tot;
      }
      ++dface_ind;
    }
    if(tot > sizeThis){ 
      numBdryFaces[elem] = tot;
      // After copying is done
      updateAdjElemFlags(dual_elems, dual_faces, elem, bdryFacesUpdate, nbf);
    }
  };

  o::LO total = 1;

  while(total){
    o::parallel_for(msh.nelems(), fill, "FillBdryFacesInADepth");

    total = o::parallel_reduce(msh.nelems(), OMEGA_H_LAMBDA(
      const int i, double& update){
        update += bdryFacesUpdate(i);  //[]
    }, "BdryFillCheck");
  }
}







} //namespace

#endif //define

