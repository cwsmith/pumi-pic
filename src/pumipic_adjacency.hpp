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


//In particles.cpp
//prevPt, finalPt, collPt, Pid, elemid,flag
//typedef MemberTypes<Vector3d, Vector3d, Vector3d, int, int, int > Particle;

OMEGA_H_INLINE bool search_mesh( SellCSigma<Particle>* scs,
  const MeshData &meshData, o::LO &loops, o::LO limit=0)
{
  const auto &side_is_exposed = meshData.side_is_exposed;
  const auto &mesh2verts = meshData.mesh2verts; 
  const auto &face_verts = meshData.face_verts; 
  const auto nelems = meshData.mesh.nelems(); //scs->num_elems ?
  const auto &down_r2fs = meshData.down_r2f.ab2b;
  const auto &dual_faces = meshData.dual.ab2b;
  const auto &dual_elems = meshData.dual.a2ab;

  const int totNumPtcls = elem_ids.size();
  o::Write<o::LO> elem_ids_next(totNumPtcls,-1);

  //particle search: adjacency + boundary crossing
  auto search_ptcl = OMEGA_H_LAMBDA( o::LO ielem)
  {
    // NOTE ielem is taken as sequential from 0 ... is it elementID ? TODO verify it
    const auto tetv2v = o::gather_verts<4>(mesh2verts, ielem);
    const auto M = o::gather_vectors<4, 3>(coords, tetv2v);

    // parallel_for loop for groups of remaining particles in this element
    //......

    // Each group of particles inside the parallel_for.
    // TODO Change ntpcl, ip start and limit. Update global(?) indices inside.
    for(o::LO ip = 0; ip < totNumPtcls; ++ip) //HACK - each element checks all particles
    {
      //skip if the particle is not in this element or has been found
      if(elem_ids[ip] != ielem || part_flags[ip] <= 0) continue;

      if(verbose > 1)
        std::cerr << "Elem " << ielem << " ptcl:" << ip << "\n";
        
      const o::Vector<3> orig{x0[ip], y0[ip], z0[ip]};
      const o::Vector<3> dest{x[ip], y[ip], z[ip]};
      
      o::Write<o::Real> bcc(4, -1.0);

      //TESTING. Check particle origin containment in current element
      find_barycentric_tet(M, orig, bcc);
      if(verbose > 2 && !(all_positive(bcc.data(), 4)))
          std::cerr << "ORIGIN ********NOT in elemet_id " << ielem << " \n";

      find_barycentric_tet(M, dest, bcc);

      //check if the destination is in this element
      if(all_positive(bcc.data(), 4, 0)) {//SURFACE_EXCLUDE)) TODO
        // TODO interpolate Fields to ptcl position, and store them, for push
        // interpolateFields(bcc, ptcls);
        elem_ids_next[ip] = elem_ids[ip];
        part_flags.data()[ip] = -1;
        if(verbose > 1) {
            std::cerr << "********found in " << ielem << " \n";
            print_matrix(M);
        }
        continue;
      }
       //get element ID
      //TODO get map from omega methods. //2,3 nodes of faces. 0,2,1; 0,1,3; 1,2,3; 2,0,3
      o::LOs fmap{2,1,1,3,2,3,0,3}; 
      auto dface_ind = (*dual_elems)[ielem];
      const auto beg_face = ielem *4;
      const auto end_face = beg_face +4;
      o::LO f_index = 0;
      bool inverse;

      for(auto iface = beg_face; iface < end_face; ++iface) {//not 0..3
        const auto face_id = down_r2fs[iface];
        if(verbose > 2)  
          std::cout << " \nFace: " << face_id << " dface_ind " <<  dface_ind << "\n";

        o::Vector<3> xpoint{0,0,0};
        auto fv2v = o::gather_verts<3>(face_verts, face_id); //Few<LO, 3>

        const auto face = o::gather_vectors<3, 3>(coords, fv2v);
        o::LO matInd1 = fmap[f_index*2], matInd2 = fmap[f_index*2+1];

        if(verbose > 2) {
          std::cout << "Face_local_index "<< fv2v.data()[0] << " " << fv2v.data()[1] << " " << fv2v.data()[2] << "\n";
          std::cout << "Mat index "<< tetv2v[matInd1] << " " << tetv2v[matInd2] << " " <<  matInd1 << " " << matInd2 << " \n";
          std::cout << "Mat dat ind " <<  tetv2v.data()[0] << " " << tetv2v.data()[1] << " "
                   << tetv2v.data()[2] << " " << tetv2v.data()[3] << "\n";
        }


        if(fv2v.data()[1] == tetv2v[matInd1] && fv2v.data()[2] == tetv2v[matInd2])
          inverse = false;
        else {// if(fv2v.data()[1] == tetv2v[matInd2] && fv2v.data()[2] == tetv2v[matInd1])
          inverse = true;
        }

        //TODO not useful
        auto fcoords = o::gather_vectors<3, 3>(coords, fv2v);
        auto base = o::simplex_basis<3, 2>(fcoords); //edgres = Matrix<2,3>
        auto snormal = o::normalize(o::cross(base[0], base[1]));

        o::LO dummy = -1;
        bool detected = line_triangle_intx_simple(face, orig, dest, xpoint, dummy, inverse);
        if(detected && verbose >1)
            std::cout << " Detected: For el=" << ielem << "\n";

        if(detected && side_is_exposed[face_id]) {
           part_flags.data()[ip] = -1;
           for(o::LO i=0; i<3; ++i)xpoints[ip*3+i] = xpoint.data()[i];
           //store current face_id and element_ids

           if(verbose > 1)
             print_osh_vector(xpoint, "COLLISION POINT");

           elem_ids_next[ip] = -1;
           break;
         }
         else if(detected && !side_is_exposed[face_id]) {
          //OMEGA_H_CHECK(side2side_elems[side + 1] - side2side_elems[side] == 2);
           auto adj_elem  = (*dual_faces)[dface_ind];
           if(verbose > 1)
             std::cout << "Deletected For el=" << ielem << " ;face_id=" << down_r2fs[iface]
                     << " ;ADJ elem= " << adj_elem << "\n";

           elem_ids_next[ip] = adj_elem;
           break;
         }

         if(!side_is_exposed[face_id]) {//TODO for DEBUG
           if(verbose > 2)
             std::cout << "adj_element_across_this_face " << (*dual_faces)[dface_ind] << "\n";
           const o::LO min_ind = min_index(bcc.data(), 4);
           if(f_index == min_ind) {
             if(verbose > 2)
               std::cout << "Min_bcc el|face_id=" << ielem << "," << down_r2fs[iface]
                     << " :unused adj_elem= " << (*dual_faces)[dface_ind] << "\n";
            if(!detected) {
              elem_ids_next[ip] = (*dual_faces)[dface_ind];
              if(verbose > 2)
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
  while(!found) {
    if(verbose > 0) 
      fprintf(stderr, "------------ %d ------------\n", loops);
    //TODO check if particle is on boundary and remove from list if so.

    // Searching all elements. TODO exclude those done ?
    o::parallel_for(nelems,  search_ptcl, "search_ptcl");
    found = true;
    auto cp_elm_ids = OMEGA_H_LAMBDA( o::LO i) {
      elem_ids[i] = elem_ids_next[i];
    };
    o::parallel_for(elem_ids.size(), cp_elm_ids, "copy_elem_ids");

    // TODO synchronize

    //TODO this could be a sequential bottle-neck
    for(int i=0; i<totNumPtcls; ++i) { 
      if(part_flags[i] > 0) {
        found = false; 
        break;
      } 
    }
    //Copy particle data from previous to next (adjacent) element
    ++loops;

    if(limit && loops > limit) 
      break;
  }

  if(verbose > 0)
    std::cerr << "search iterations " << loops << "\n";

  return found;
} //search_mesh

} //namespace

#endif //define

