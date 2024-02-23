use std::fs::File;

use num_traits::AsPrimitive;

use crate::topology::{DynamicTriangle, DynamicVertex};

// --------------------------------------------

#[allow(clippy::identity_op)]
pub fn make_super_triangle<T>(
    vtx2xy: &mut Vec<nalgebra::Vector2<T>>,
    vmin: &[T; 2],
    vmax: &[T; 2]) -> (Vec<DynamicTriangle>, Vec<DynamicVertex>)
    where T: num_traits::Float + 'static + std::fmt::Debug + std::default::Default,
          f64: AsPrimitive<T>
{ // super triangle
    let mut vtx2tri = vec!(DynamicVertex { e: usize::MAX, d: 0 }; vtx2xy.len());
    //
    assert_eq!(vtx2tri.len(), vtx2xy.len());
    let (max_len, center) = {
        let vsize = [vmax[0] - vmin[0], vmax[1] - vmin[1]];
        let max_len = if vsize[0] > vsize[1] { vsize[0] } else { vsize[1] };
        (max_len, [(vmin[0] + vmax[0]) * 0.5_f64.as_(), (vmin[1] + vmax[1]) * 0.5_f64.as_()])
    };
    let tri_len: T = max_len * 4_f64.as_();
    let tmp_len: T = tri_len * (3.0_f64.sqrt() / 6.0_f64).as_();
    let npo = vtx2xy.len();
    //
    vtx2xy.resize(npo + 3, Default::default());
    vtx2xy[npo + 0] = nalgebra::Vector2::<T>::new(center[0], center[1] + 2_f64.as_() * tmp_len);
    vtx2xy[npo + 1] = nalgebra::Vector2::<T>::new(center[0] - 0.5_f64.as_() * tri_len, center[1] - tmp_len);
    vtx2xy[npo + 2] = nalgebra::Vector2::<T>::new(center[0] + 0.5_f64.as_() * tri_len, center[1] - tmp_len);
    //
    vtx2tri.resize(npo + 3, DynamicVertex { e: 0, d: 0 });
    vtx2tri[npo + 0].e = 0;
    vtx2tri[npo + 0].d = 0;
    vtx2tri[npo + 1].e = 0;
    vtx2tri[npo + 1].d = 1;
    vtx2tri[npo + 2].e = 0;
    vtx2tri[npo + 2].d = 2;
    //
    let mut tri2vtx = vec!(DynamicTriangle { v: [0; 3], s: [0; 3] }; 1);
    {
        let tri = &mut tri2vtx[0];
        tri.v = [npo + 0, npo + 1, npo + 2];
        tri.s = [usize::MAX; 3];
    }
    (tri2vtx, vtx2tri)
}


pub fn add_points_to_mesh<T>(
    tri2vtx: &mut Vec<DynamicTriangle>,
    vtx2tri: &mut Vec<DynamicVertex>,
    vtx2xy: &Vec<nalgebra::Vector2<T>>,
    i_vtx: usize)
    where T: num_traits::Float + 'static + Copy + std::fmt::Debug,
          f64: num_traits::AsPrimitive<T>
{
    assert_eq!(vtx2xy.len(), vtx2tri.len());
    if vtx2tri[i_vtx].e != usize::MAX { return; } // already added
    let po_add = vtx2xy[i_vtx];
    for itri in 0..tri2vtx.len() {
        let areas = [
            del_geo::tri2::area(&po_add, &vtx2xy[tri2vtx[itri].v[1]], &vtx2xy[tri2vtx[itri].v[2]]),
            del_geo::tri2::area(&po_add, &vtx2xy[tri2vtx[itri].v[2]], &vtx2xy[tri2vtx[itri].v[0]]),
            del_geo::tri2::area(&po_add, &vtx2xy[tri2vtx[itri].v[0]], &vtx2xy[tri2vtx[itri].v[1]])];
        let area_sum: T = areas[0] + areas[1] + areas[2];
        assert!(area_sum > T::zero());
        let (&area_min, iedge) = areas.iter().zip(0..)
            .min_by(|a, b| a.0.partial_cmp(b.0).expect("NaN area"))
            .unwrap();
        // dbg!(area_sum, areas, area_min, iedge);
        if area_min <= -area_sum * 1.0e-10f64.as_() { continue; } // the point is out of the triangle
        //
        if area_min > area_sum * 1.0e-3f64.as_() {
            crate::topology::insert_a_point_inside_an_element(i_vtx, itri, vtx2tri, tri2vtx);
        } else {
            crate::topology::insert_point_on_elem_edge(i_vtx, itri, iedge, vtx2tri, tri2vtx);
        }
        return;
    }
    panic!();
}

pub fn should_flip<T>(
    i_tri0: usize,
    i_node0: usize,
    tri2vtx: &Vec<DynamicTriangle>,
    vtx2xy: &[nalgebra::Vector2<T>]) -> bool
    where T: num_traits::Float + std::fmt::Debug + 'static + Copy,
          f64: AsPrimitive<T>
{
    if tri2vtx[i_tri0].s[i_node0] >= tri2vtx.len() { return false; }// there is adjacent triangle
    let j_tri0 = tri2vtx[i_tri0].s[i_node0];
    let j_node0 = crate::topology::find_adjacent_edge_index(
        &tri2vtx[i_tri0], i_node0, tri2vtx);
    assert_eq!(tri2vtx[j_tri0].s[j_node0], i_tri0);
    let pj0 = vtx2xy[tri2vtx[j_tri0].v[j_node0]];
    let pi0 = vtx2xy[tri2vtx[i_tri0].v[i_node0]];
    let pi1 = vtx2xy[tri2vtx[i_tri0].v[(i_node0 + 1) % 3]];
    let pi2 = vtx2xy[tri2vtx[i_tri0].v[(i_node0 + 2) % 3]];
    let a_i0_i1_i2 = del_geo::tri2::area(&pi0, &pi1, &pi2);
    let a_j0_i2_i1 = del_geo::tri2::area(&pj0, &pi2, &pi1);
    assert!(a_i0_i1_i2 > T::zero());
    assert!(a_j0_i2_i1 > T::zero());
    let area_diamond = a_i0_i1_i2 + a_j0_i2_i1;
    let a_i0_i1_j0 = del_geo::tri2::area(&pi0, &pi1, &pj0);
    let a_i0_j0_i2 = del_geo::tri2::area(&pi0, &pj0, &pi2);
    if a_i0_i1_j0 < area_diamond * 1.0e-3f64.as_() { return false; }
    if a_i0_j0_i2 < area_diamond * 1.0e-3f64.as_() { return false; }
    let cc = del_geo::tri2::circumcenter(&pi0, &pi1, &pi2);
    let rad = del_geo::edge2::length_squared(&cc, &pi0);
    let dist = del_geo::edge2::length_squared(&cc, &pj0);
    if dist >= rad { return false; }
    true
}

pub fn delaunay_around_point<T>(
    i_vtx0: usize,
    vtx2tri: &mut Vec<DynamicVertex>,
    tri2vtx: &mut Vec<DynamicTriangle>,
    vtx2xy: &Vec<nalgebra::Vector2<T>>)
    where T: num_traits::Float + 'static + Copy + std::fmt::Debug,
          f64: AsPrimitive<T>
{
    use crate::topology::{flip_edge, move_ccw, move_cw};
    assert_eq!(vtx2xy.len(), vtx2tri.len());
    assert!(i_vtx0 < vtx2tri.len());
    if vtx2tri[i_vtx0].e == usize::MAX { return; }

    let mut i_tri0 = vtx2tri[i_vtx0].e;
    let mut i_node0 = vtx2tri[i_vtx0].d;

    // ---------------------------
    // go counter-clock-wise
    let mut flag_is_wall = false;
    loop {
        assert_eq!(tri2vtx[i_tri0].v[i_node0], i_vtx0);
        if should_flip(i_tri0, i_node0, tri2vtx, vtx2xy) { // there is adjacent triangle
            flip_edge(i_tri0, i_node0, vtx2tri, tri2vtx); // this edge is not on the edge and should be successful
            i_node0 = 2;
            assert_eq!(tri2vtx[i_tri0].v[i_node0], i_vtx0); // this is the rule from FlipEdge function
            continue; // need to check the fliped element
        }
        if !move_ccw(&mut i_tri0, &mut i_node0, usize::MAX, tri2vtx) {
            flag_is_wall = true;
            break;
        }
        if i_tri0 == vtx2tri[i_vtx0].e {
            break;
        }
    }
    if !flag_is_wall { return; }

    // ----------------------------
    // go clock-wise
    loop {
        assert_eq!(tri2vtx[i_tri0].v[i_node0], i_vtx0);
        if should_flip(i_tri0, i_node0, tri2vtx, vtx2xy) {
            let j_tri0 = tri2vtx[i_tri0].s[i_node0];
            flip_edge(i_tri0, i_node0, vtx2tri, tri2vtx);
            i_tri0 = j_tri0;
            i_node0 = 1;
            assert_eq!(tri2vtx[i_tri0].v[i_node0], i_vtx0);
            continue;
        }
        if !move_cw(&mut i_tri0, &mut i_node0, usize::MAX, tri2vtx) { return; }
    }
}

pub fn meshing_initialize(
    vtx2xy: &mut Vec<nalgebra::Vector2<f32>>) -> (Vec<DynamicTriangle>, Vec<DynamicVertex>)
{
    let (mut tri2vtx, mut vtx2tri) = {
        let aabb = del_geo::aabb2::from_vtx2vec(vtx2xy);
        make_super_triangle(
            vtx2xy,
            &[aabb[0], aabb[1]], &[aabb[2], aabb[3]])
    };
    for i_vtx in 0..vtx2tri.len() - 3 {
        add_points_to_mesh(
            &mut tri2vtx, &mut vtx2tri, vtx2xy,
            i_vtx);
        // crate::trimesh2_util::save_tri_mesh(format!("vtx_{}a.obj",ip), &tris, vtx2xy);
        delaunay_around_point(
            i_vtx,
            &mut vtx2tri, &mut tri2vtx, vtx2xy);
        // crate::trimesh2_util::save_tri_mesh(format!("vtx_{}b.obj",ip), &tris, vtx2xy);
    }
    (tri2vtx, vtx2tri)
}


fn find_edge_point_across_edge<T>(
    itri0: &mut usize,
    inotri0: &mut usize,
    inotri1: &mut usize,
    ratio: &mut T,
    ipo0: usize,
    ipo1: usize,
    vtx_tri: &[DynamicVertex],
    tri_vtx: &[DynamicTriangle],
    vtx_xy: &[nalgebra::Vector2<T>]) -> bool
    where T: num_traits::Float + 'static + Copy + std::fmt::Debug,
          f64: AsPrimitive<T>
{
    use crate::topology::find_adjacent_edge_index;
    let itri_ini = vtx_tri[ipo0].e;
    let inotri_ini = vtx_tri[ipo0].d;
    let mut inotri_cur = inotri_ini;
    let mut itri_cur = itri_ini;
    loop {
        assert_eq!(tri_vtx[itri_cur].v[inotri_cur], ipo0);
        {
            let inotri2 = (inotri_cur + 1) % 3;
            let inotri3 = (inotri_cur + 2) % 3;
            let area0 = del_geo::tri2::area(
                &vtx_xy[ipo0],
                &vtx_xy[tri_vtx[itri_cur].v[inotri2]],
                &vtx_xy[ipo1]);
            if area0 > -(1.0e-20_f64.as_()) {
                let area1 = del_geo::tri2::area(
                    &vtx_xy[ipo0],
                    &vtx_xy[ipo1],
                    &vtx_xy[tri_vtx[itri_cur].v[inotri3]]);
                if area1 > -(1.0e-20_f64.as_()) {
                    dbg!(area0,area1);
                    assert!(area0 + area1 > 1.0e-20_f64.as_());
                    *ratio = area0 / (area0 + area1);
                    *itri0 = itri_cur;
                    *inotri0 = inotri2;
                    *inotri1 = inotri3;
                    return true;
                }
            }
        }
        {
            let inotri2 = (inotri_cur + 1) % 3;
            let itri_nex = tri_vtx[itri_cur].s[inotri2];
            if itri_nex == usize::MAX { break; }
            let jnob = find_adjacent_edge_index(
                &tri_vtx[itri_nex], inotri2, tri_vtx);
            let inotri3 = (jnob + 1) % 3;
            assert!(itri_nex < tri_vtx.len());
            assert_eq!(tri_vtx[itri_nex].v[inotri3], ipo0);
            if itri_nex == itri_ini {
                *itri0 = 0;
                *inotri0 = 0;
                *inotri1 = 0;
                *ratio = 0_f64.as_();
                return false;
            }
            itri_cur = itri_nex;
            inotri_cur = inotri3;
        }
    }

    inotri_cur = inotri_ini;
    itri_cur = itri_ini;
    loop {
        assert_eq!(tri_vtx[itri_cur].v[inotri_cur], ipo0);
        {
            let inotri2 = (inotri_cur + 1) % 3; // indexRot3[1][inotri_cur];
            let inotri3 = (inotri_cur + 2) % 3; // indexRot3[2][inotri_cur];
            let area0 = del_geo::tri2::area(
                &vtx_xy[ipo0],
                &vtx_xy[tri_vtx[itri_cur].v[inotri2]],
                &vtx_xy[ipo1]);
            if area0 > -(1.0e-20_f64.as_()) {
                let area1 = del_geo::tri2::area(
                    &vtx_xy[ipo0],
                    &vtx_xy[ipo1],
                    &vtx_xy[tri_vtx[itri_cur].v[inotri3]]);
                if area1 > -(1.0e-20_f64.as_()) {
                    assert!(area0 + area1 > 1.0e-20_f64.as_());
                    *ratio = area0 / (area0 + area1);
                    *itri0 = itri_cur;
                    *inotri0 = inotri2;
                    *inotri1 = inotri3;
                    return true;
                }
            }
        }
        {
            let inotri2 = (inotri_cur + 2) % 3;
            let itri_nex = tri_vtx[itri_cur].s[inotri2];
            let jnob = find_adjacent_edge_index(&tri_vtx[itri_cur], inotri2, tri_vtx);
            let inotri3 = (jnob + 1) % 3;
            assert_eq!(tri_vtx[itri_nex].v[inotri3], ipo0);
            if itri_nex == itri_ini {
                panic!();
            }
            itri_cur = itri_nex;
            inotri_cur = inotri3;
        }
    }
}

pub fn enforce_edge<T>(
    vtx2tri: &mut Vec<DynamicVertex>,
    tri2vtx: &mut Vec<DynamicTriangle>,
    i0_vtx: usize,
    i1_vtx: usize,
    vtx2xy: &[nalgebra::Vector2<T>])
    where T: num_traits::Float + 'static + Copy + std::fmt::Debug + std::fmt::Display,
          f64: AsPrimitive<T>
{
    assert_eq!(vtx2xy.len(), vtx2tri.len());
    use crate::topology::{
        find_adjacent_edge_index,
        find_edge_by_looking_around_point,
        flip_edge};
    assert!(i0_vtx < vtx2tri.len());
    assert!(i1_vtx < vtx2tri.len());
    loop {
        let mut itri0: usize = usize::MAX;
        let mut inotri0: usize = 0;
        let mut inotri1: usize = 0;
        if find_edge_by_looking_around_point(
            &mut itri0, &mut inotri0, &mut inotri1,
            i0_vtx, i1_vtx,
            vtx2tri, tri2vtx) { // this edge divides outside and inside
            assert_ne!(inotri0, inotri1);
            assert!(inotri0 < 3);
            assert!(inotri1 < 3);
            assert_eq!(tri2vtx[itri0].v[inotri0], i0_vtx);
            assert_eq!(tri2vtx[itri0].v[inotri1], i1_vtx);
            let ied0 = 3 - inotri0 - inotri1;
            {
                let itri1 = tri2vtx[itri0].s[ied0];
                let ied1 = find_adjacent_edge_index(&tri2vtx[itri0], ied0, tri2vtx);
                assert_eq!(tri2vtx[itri1].s[ied1], itri0);
                tri2vtx[itri1].s[ied1] = usize::MAX;
                tri2vtx[itri0].s[ied0] = usize::MAX;
            }
            break;
        } else { // this edge is divided from connection outer triangle
            let mut ratio: T = 0_f64.as_();
            if !find_edge_point_across_edge(
                &mut itri0, &mut inotri0, &mut inotri1, &mut ratio,
                i0_vtx, i1_vtx,
                vtx2tri, tri2vtx, vtx2xy) { panic!(); }
            assert!(ratio > -(1.0e-20_f64.as_()));
            // assert!( ratio < 1_f64.as_() + 1.0e-20_f64.as_());
            assert!(del_geo::tri2::area(&vtx2xy[i0_vtx], &vtx2xy[tri2vtx[itri0].v[inotri0]], &vtx2xy[i1_vtx]) > 1.0e-20_f64.as_());
            assert!(del_geo::tri2::area(&vtx2xy[i0_vtx], &vtx2xy[i1_vtx], &vtx2xy[tri2vtx[itri0].v[inotri1]]) > 1.0e-20_f64.as_());
//            std::cout << ratio << std::endl;
            if ratio < 1.0e-20_f64.as_() {
                panic!();
            } else if ratio > 1.0_f64.as_() - 1.0e-10_f64.as_() {
                panic!();
            } else {
                let ied0 = 3 - inotri0 - inotri1;
                assert!(tri2vtx[itri0].s[ied0] < tri2vtx.len());
                let res = flip_edge(itri0, ied0, vtx2tri, tri2vtx);
                if !res {
                    break;
                }
            }
        }
    }
}


pub fn delete_unreferenced_points(
    vtx_xy: &mut Vec<nalgebra::Vector2<f32>>,
    vtx_tri: &mut Vec<DynamicVertex>,
    tri_vtx: &mut [DynamicTriangle],
    point_idxs_to_delete: &Vec<usize>) {
    assert_eq!(vtx_tri.len(), vtx_xy.len());
    let mut map_po_del = Vec::<usize>::new();
    let mut npo_pos;
    {
        map_po_del.resize(vtx_tri.len(), usize::MAX - 1);
        for ipo in point_idxs_to_delete {
            map_po_del[*ipo] = usize::MAX;
        }
        npo_pos = 0;
        for ipo in 0..vtx_tri.len() {
            if map_po_del[ipo] == usize::MAX {
                continue;
            }
            map_po_del[ipo] = npo_pos;
            npo_pos += 1;
        }
    }
    {
        let vtx_tri_tmp = vtx_tri.clone();
        let vtx_xy_tmp = vtx_xy.clone();
        vtx_tri.resize(npo_pos, DynamicVertex { e: 0, d: 0 });
        vtx_xy.resize(npo_pos, Default::default());
        for ipo in 0..map_po_del.len() {
            if map_po_del[ipo] == usize::MAX {
                continue;
            }
            let ipo1 = map_po_del[ipo];
            vtx_tri[ipo1] = vtx_tri_tmp[ipo].clone();
            vtx_xy[ipo1] = vtx_xy_tmp[ipo];
        }
    }
    for (itri, tri) in tri_vtx.iter_mut().enumerate() {
        for ifatri in 0..3 {
            let ipo = tri.v[ifatri];
            assert_ne!(map_po_del[ipo], usize::MAX);
            tri.v[ifatri] = map_po_del[ipo];
            vtx_tri[ipo].e = itri;
            vtx_tri[ipo].d = ifatri;
        }
    }
}


pub fn laplacian_mesh_smoothing_around_point<T>(
    vtx2xy: &mut Vec<nalgebra::Vector2<T>>,
    i_vtx0: usize,
    vtx2tri: &Vec<DynamicVertex>,
    tri2vtx: &Vec<DynamicTriangle>) -> bool
    where T: num_traits::Float + Clone + std::fmt::Debug + std::ops::AddAssign + std::ops::DivAssign + 'static + Copy,
          usize: num_traits::AsPrimitive<T>,
          f64: AsPrimitive<T>
{
    use crate::topology::move_ccw;
    assert_eq!(vtx2xy.len(), vtx2tri.len());
    let mut i_tri0 = vtx2tri[i_vtx0].e;
    let mut i_node0 = vtx2tri[i_vtx0].d;
    let pos_before = vtx2xy[i_vtx0];
    let mut pos_new = vtx2xy[i_vtx0];
    let mut ntri_around: usize = 1;
    loop { // counter-clock wise
        assert!(i_tri0 < tri2vtx.len() && i_node0 < 3 && tri2vtx[i_tri0].v[i_node0] == i_vtx0);
        pos_new += vtx2xy[tri2vtx[i_tri0].v[(i_node0 + 1) % 3]];
        ntri_around += 1;
        if !move_ccw(&mut i_tri0, &mut i_node0, usize::MAX, tri2vtx) { return false; }
        if i_tri0 == vtx2tri[i_vtx0].e { break; }
    }
    vtx2xy[i_vtx0] = pos_new / ntri_around.as_();
    //
    let mut flipped = false;
    i_tri0 = vtx2tri[i_vtx0].e;
    i_node0 = vtx2tri[i_vtx0].d;
    loop { // counter-clock wise
        let area = area_of_triangle_in_mesh(tri2vtx, vtx2xy, i_tri0);
        if area < T::zero() {
            flipped = true;
            break;
        }
        assert!(i_tri0 < tri2vtx.len() && i_node0 < 3 && tri2vtx[i_tri0].v[i_node0] == i_vtx0);
        if !move_ccw(&mut i_tri0, &mut i_node0, usize::MAX, tri2vtx) { return false; }
        if i_tri0 == vtx2tri[i_vtx0].e { break; }
    }
    if flipped {
        vtx2xy[i_vtx0] = pos_before;
    }
    true
}

#[allow(clippy::identity_op)]
pub fn save_tri_mesh<P, T>(
    filepath: P,
    tri2vtx: &[DynamicTriangle],
    vtx2xy: &[nalgebra::Vector2<T>])
    where P: AsRef<std::path::Path>,
          T: std::fmt::Display
{
    use std::io::Write;
    let file = File::create(filepath).expect("file not found.");
    let mut file = std::io::BufWriter::new(file);
    for vtx in vtx2xy {
        writeln!(file, "v {} {} {}",
                 vtx[0], vtx[1], 0.0).expect("fail");
    }
    for tri in tri2vtx {
        writeln!(file, "f {} {} {}",
                 tri.v[0] + 1,
                 tri.v[1] + 1,
                 tri.v[2] + 1).expect("fail");
    }
}

fn area_of_triangle_in_mesh<T>(
    tri2vtx: &[DynamicTriangle],
    vtx2xy: &[nalgebra::Vector2<T>],
    i_tri: usize) -> T
    where T: num_traits::Float + 'static + Copy,
          f64: AsPrimitive<T>
{
    let i0 = tri2vtx[i_tri].v[0];
    let i1 = tri2vtx[i_tri].v[1];
    let i2 = tri2vtx[i_tri].v[2];
    del_geo::tri2::area(&vtx2xy[i0], &vtx2xy[i1], &vtx2xy[i2])
}