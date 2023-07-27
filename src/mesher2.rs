//! 2D dynamic mesh editing utility code

use num_traits::AsPrimitive;
use crate::topology::{DynamicTriangle, DynamicVertex};

fn bounding_box2<VEC>(
    vtx_xy: &mut Vec<VEC>)
    -> (VEC, VEC)
    where VEC: Copy + std::ops::IndexMut<usize>,
          <VEC as std::ops::Index<usize>>::Output: PartialOrd + Sized + Copy
{
    let (mut vmin, mut vmax) = (vtx_xy[0], vtx_xy[0]);
    for ivtx in 1..vtx_xy.len() {
        let xy = &vtx_xy[ivtx];
        if xy[0] < vmin[0] { vmin[0] = xy[0]; }
        if xy[0] > vmax[0] { vmax[0] = xy[0]; }
        if xy[1] < vmin[1] { vmin[1] = xy[1]; }
        if xy[1] > vmax[1] { vmax[1] = xy[1]; }
    }
    (vmin, vmax)
}

fn area_tri2<T>(
    v1: &nalgebra::Vector2<T>,
    v2: &nalgebra::Vector2<T>,
    v3: &nalgebra::Vector2<T>) -> T
    where T: num_traits::Float + 'static + Copy,
          f64: num_traits::AsPrimitive<T>
{
    ((v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1])) / 2_f64.as_()
}

fn squared_distance<T>(
    v1: &nalgebra::Vector2<T>,
    v2: &nalgebra::Vector2<T>) -> T
    where T: num_traits::real::Real + 'static + Copy,
{
    (v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[1] - v2[1]) * (v1[1] - v2[1])
}

fn det_delaunay<T>(
    p0: &nalgebra::Vector2<T>,
    p1: &nalgebra::Vector2<T>,
    p2: &nalgebra::Vector2<T>,
    p3: &nalgebra::Vector2<T>) -> i32
    where T: num_traits::Float + 'static + Copy,
          f64: num_traits::AsPrimitive<T>
{
    let area = area_tri2(p0, p1, p2);
    if area.abs() < 1.0e-10_f64.as_() {
        return 3;
    }
    let tmp_val = 1_f64.as_() / (area * area * 16_f64.as_());

    let dtmp0 = squared_distance(p1, p2);
    let dtmp1 = squared_distance(p0, p2);
    let dtmp2 = squared_distance(p0, p1);

    let etmp0: T = tmp_val * dtmp0 * (dtmp1 + dtmp2 - dtmp0);
    let etmp1: T = tmp_val * dtmp1 * (dtmp0 + dtmp2 - dtmp1);
    let etmp2: T = tmp_val * dtmp2 * (dtmp0 + dtmp1 - dtmp2);

    let out_center = nalgebra::Vector2::<T>::new(
        etmp0 * p0[0] + etmp1 * p1[0] + etmp2 * p2[0],
        etmp0 * p0[1] + etmp1 * p1[1] + etmp2 * p2[1]);

    let qradius = squared_distance(&out_center, &p0);
    let qdistance = squared_distance(&out_center, &p3);

//	assert( fabs( qradius - SquareLength(out_center,p1) ) < 1.0e-10*qradius );
//	assert( fabs( qradius - SquareLength(out_center,p2) ) < 1.0e-10*qradius );

    let tol: T = 1.0e-20.as_();
    return if qdistance > qradius * (1_f64.as_() + tol) { 2 }    // outside the circumcircle
    else {
        if qdistance < qradius * (1_f64.as_() - tol) { 0 }    // inside the circumcircle
        else { 1 }    // on the circumcircle
    };
}

// --------------------------------------------

pub fn make_super_triangle<T>(
    tri2vtx: &mut Vec<DynamicTriangle>,
    vtx2tri: &mut Vec<DynamicVertex>,
    vtx2xy: &mut Vec<nalgebra::Vector2<T>>,
    vmin: &[T; 2],
    vmax: &[T; 2])
    where T: num_traits::Float + 'static + std::fmt::Debug + std::default::Default,
          f64: AsPrimitive<T>
{ // super triangle
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
    tri2vtx.clear();
    tri2vtx.resize(1, DynamicTriangle { v: [0; 3], s: [0; 3] });
    let tri = &mut tri2vtx[0];
    tri.v = [npo + 0, npo + 1, npo + 2];
    tri.s = [usize::MAX; 3];
}


pub fn add_points_to_mesh<T>(
    tris: &mut Vec<DynamicTriangle>,
    vtx2tri: &mut Vec<DynamicVertex>,
    vtx2xy: &Vec<nalgebra::Vector2<T>>,
    idx_vtx: usize,
    min_tri_area: T)
    where T: num_traits::Float + 'static + Copy,
          f64: num_traits::AsPrimitive<T>
{
    use crate::topology::{
        find_adjacent_edge_index,
        insert_a_point_inside_an_element,
        insert_point_on_elem_edge,
    };

    assert_eq!(vtx2xy.len(), vtx2tri.len());
    if vtx2tri[idx_vtx].e != usize::MAX { return; } // already added
    let po_add = vtx2xy[idx_vtx];
    let mut itri_in = usize::MAX;
    let mut iedge = usize::MAX;
    let mut iflg1;
    let mut iflg2;
    for itri in 0..tris.len() {
        iflg1 = 0;
        iflg2 = 0;
        let a0 = area_tri2(&po_add, &vtx2xy[tris[itri].v[1]], &vtx2xy[tris[itri].v[2]]);
        let a1 = area_tri2(&po_add, &vtx2xy[tris[itri].v[2]], &vtx2xy[tris[itri].v[0]]);
        let a2 = area_tri2(&po_add, &vtx2xy[tris[itri].v[0]], &vtx2xy[tris[itri].v[1]]);
        if a0 > min_tri_area {
            iflg1 += 1;
            iflg2 += 0;
        }
        if a1 > min_tri_area {
            iflg1 += 1;
            iflg2 += 1;
        }
        if a2 > min_tri_area {
            iflg1 += 1;
            iflg2 += 2;
        }
        if iflg1 == 3 { // add in triangle
            itri_in = itri;
            break;
        } else if iflg1 == 2 {
            // add in edge
            let ied0 = 3 - iflg2;
            let ipo_e0 = tris[itri].v[(ied0 + 1) % 3];
            let ipo_e1 = tris[itri].v[(ied0 + 2) % 3];
            let itri_s = tris[itri].s[ied0];
            if itri_s == usize::MAX { return; }
            let jno0 = find_adjacent_edge_index(&tris[itri], ied0, tris);
            assert_eq!(tris[itri_s].v[(jno0 + 2) % 3], ipo_e0);
            assert_eq!(tris[itri_s].v[(jno0 + 1) % 3], ipo_e1);
            let inoel_d = jno0;
            assert_eq!(tris[itri_s].s[inoel_d], itri);
            let ipo_d = tris[itri_s].v[inoel_d];
            assert!(area_tri2(&po_add, &vtx2xy[ipo_e1], &vtx2xy[tris[itri].v[ied0]]) > min_tri_area);
            assert!(area_tri2(&po_add, &vtx2xy[tris[itri].v[ied0]], &vtx2xy[ipo_e0]) > min_tri_area);
            if area_tri2(&po_add, &vtx2xy[ipo_e0], &vtx2xy[ipo_d]) < min_tri_area { continue; }
            if area_tri2(&po_add, &vtx2xy[ipo_d], &vtx2xy[ipo_e1]) < min_tri_area { continue; }
            let det_d = det_delaunay(
                &po_add,
                &vtx2xy[ipo_e0], &vtx2xy[ipo_e1], &vtx2xy[ipo_d]);
            if det_d == 2 || det_d == 1 {
                continue;
            }
            itri_in = itri;
            iedge = ied0;
            break;
        }
    }
    if itri_in == usize::MAX {
        //std::cout << "super triangle failure " << iflg1 << " " << iflg2 << std::endl;
        panic!();
    }
    if iedge == usize::MAX {
        insert_a_point_inside_an_element(idx_vtx, itri_in, vtx2tri, tris);
    } else {
        insert_point_on_elem_edge(idx_vtx, itri_in, iedge, vtx2tri, tris);
    }
}

pub fn delaunay_around_point<T>(
    ipo0: usize,
    vtx2tri: &mut Vec<DynamicVertex>,
    tris: &mut Vec<DynamicTriangle>,
    vtx2xy: &Vec<nalgebra::Vector2<T>>)
    where T: num_traits::Float + 'static + Copy,
          f64: AsPrimitive<T>
{
    use crate::topology::{find_adjacent_edge_index, flip_edge, move_ccw, move_cw};
    assert_eq!(vtx2xy.len(), vtx2tri.len());
    assert!(ipo0 < vtx2tri.len());
    if vtx2tri[ipo0].e == usize::MAX { return; }

    let mut itri0 = vtx2tri[ipo0].e;
    let mut ino0 = vtx2tri[ipo0].d;

    // ---------------------------
    // go counter-clock-wise
    let mut flag_is_wall = false;
    loop {
        assert!(itri0 < tris.len() && ino0 < 3 && tris[itri0].v[ino0] == ipo0);
        if tris[itri0].s[ino0] < tris.len() {
            let jtri0 = tris[itri0].s[ino0];
            let jno0 = find_adjacent_edge_index(&tris[itri0], ino0, &tris);
            assert_eq!(tris[jtri0].s[jno0], itri0);
            let jpo0 = tris[jtri0].v[jno0];
            let ires = det_delaunay(
                &vtx2xy[tris[itri0].v[0]],
                &vtx2xy[tris[itri0].v[1]],
                &vtx2xy[tris[itri0].v[2]],
                &vtx2xy[jpo0]);
            if ires == 0 {
                flip_edge(itri0, ino0, vtx2tri, tris); // this edge is not on the edge and should be successfull
                ino0 = 2;
                assert_eq!(tris[itri0].v[ino0], ipo0); // this is the rule from FlipEdge function
                continue; // need to check the fliped element
            }
        }
        if !move_ccw(&mut itri0, &mut ino0, usize::MAX, tris) {
            flag_is_wall = true;
            break;
        }
        if itri0 == vtx2tri[ipo0].e {
            break;
        }
    }
    if !flag_is_wall { return; }

    // ----------------------------
    // go clock-wise

    loop {
        assert!(itri0 < tris.len() && ino0 < 3 && tris[itri0].v[ino0] == ipo0);
        if tris[itri0].s[ino0] < tris.len() {
            let jtri0 = tris[itri0].s[ino0];
            let jno0 = find_adjacent_edge_index(&tris[itri0], ino0, &tris);
            assert_eq!(tris[jtri0].s[jno0], itri0);
            let ipo_dia = tris[jtri0].v[jno0];
            let ires = det_delaunay(
                &vtx2xy[tris[itri0].v[0]],
                &vtx2xy[tris[itri0].v[1]],
                &vtx2xy[tris[itri0].v[2]],
                &vtx2xy[ipo_dia]);
            if ires == 0 { // Delaunay condition is not satisfiled
                flip_edge(itri0, ino0, vtx2tri, tris);
                itri0 = jtri0;
                ino0 = 1;
                assert_eq!(tris[itri0].v[ino0], ipo0);
                continue;
            }
        }
        if !move_cw(&mut itri0, &mut ino0, usize::MAX, tris) { return; }
    }
}

pub fn meshing_initialize(
    tri_vtx: &mut Vec<DynamicTriangle>,
    vtx_tri: &mut Vec<DynamicVertex>,
    vtx_xy: &mut Vec<nalgebra::Vector2<f32>>) {
    vtx_tri.clear();
    vtx_tri.resize(vtx_xy.len(), DynamicVertex { e: usize::MAX, d: 0 });
    {
        let (vmin, vmax) = bounding_box2::<nalgebra::Vector2<f32>>(vtx_xy);
        make_super_triangle(
            tri_vtx, vtx_tri, vtx_xy,
            &[vmin[0], vmin[1]], &[vmax[0], vmax[1]]);
    }
    {
        const MIN_TRI_AREA: f32 = 1.0e-10;
        for ip in 0..vtx_tri.len() - 3 {
            add_points_to_mesh(
                tri_vtx, vtx_tri, vtx_xy,
                ip,
                MIN_TRI_AREA);
            delaunay_around_point(
                ip,
                vtx_tri, tri_vtx, vtx_xy);
        }
    }
}


fn find_edge_point_across_edge<T>(
    itri0: &mut usize,
    inotri0: &mut usize,
    inotri1: &mut usize,
    ratio: &mut T,
    ipo0: usize,
    ipo1: usize,
    vtx_tri: &Vec<DynamicVertex>,
    tri_vtx: &Vec<DynamicTriangle>,
    vtx_xy: &Vec<nalgebra::Vector2<T>>) -> bool
    where T: num_traits::Float + 'static + Copy,
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
            let area0 = area_tri2(&vtx_xy[ipo0],
                                  &vtx_xy[tri_vtx[itri_cur].v[inotri2]],
                                  &vtx_xy[ipo1]);
            if area0 > -1.0e-20_f64.as_() {
                let area1 = area_tri2(&vtx_xy[ipo0],
                                      &vtx_xy[ipo1],
                                      &vtx_xy[tri_vtx[itri_cur].v[inotri3]]);
                if area1 > -1.0e-20_f64.as_() {
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
            let jnob = find_adjacent_edge_index(&tri_vtx[itri_nex], inotri2, tri_vtx);
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
            let area0 = area_tri2(&vtx_xy[ipo0],
                                  &vtx_xy[tri_vtx[itri_cur].v[inotri2]],
                                  &vtx_xy[ipo1]);
            if area0 > -1.0e-20_f64.as_() {
                let area1 = area_tri2(&vtx_xy[ipo0],
                                      &vtx_xy[ipo1],
                                      &vtx_xy[tri_vtx[itri_cur].v[inotri3]]);
                if area1 > -1.0e-20_f64.as_() {
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
            let jnob = find_adjacent_edge_index(&tri_vtx[itri_cur], inotri2, &tri_vtx);
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
    tris: &mut Vec<DynamicTriangle>,
    i0_vtx: usize,
    i1_vtx: usize,
    vtx2xy: &Vec<nalgebra::Vector2<T>>)
    where T: num_traits::Float + 'static + Copy,
          f64: AsPrimitive<T>
{
    use crate::topology::{
        flip_edge,
        find_edge_by_looking_around_point,
        find_adjacent_edge_index};
    assert!(i0_vtx < vtx2tri.len());
    assert!(i1_vtx < vtx2tri.len());
    loop {
        let mut itri0: usize = usize::MAX;
        let mut inotri0: usize = 0;
        let mut inotri1: usize = 0;
        if find_edge_by_looking_around_point(
            &mut itri0, &mut inotri0, &mut inotri1,
            i0_vtx, i1_vtx,
            &vtx2tri, &tris) { // this edge divide outside and inside
            assert_ne!(inotri0, inotri1);
            assert!(inotri0 < 3);
            assert!(inotri1 < 3);
            assert_eq!(tris[itri0].v[inotri0], i0_vtx);
            assert_eq!(tris[itri0].v[inotri1], i1_vtx);
            let ied0 = 3 - inotri0 - inotri1;
            {
                let itri1 = tris[itri0].s[ied0];
                let ied1 = find_adjacent_edge_index(&tris[itri0], ied0, &tris);
                assert_eq!(tris[itri1].s[ied1], itri0);
                tris[itri1].s[ied1] = usize::MAX;
                tris[itri0].s[ied0] = usize::MAX;
            }
            break;
        } else { // this edge is devided from connection outer triangle
            let mut ratio: T = 0_f64.as_();
            if !find_edge_point_across_edge(
                &mut itri0, &mut inotri0, &mut inotri1, &mut ratio,
                i0_vtx, i1_vtx,
                &vtx2tri, &tris, &vtx2xy) { panic!(); }
            assert!(ratio > -1.0e-20_f64.as_() && ratio < 1_f64.as_() + 1.0e-20_f64.as_());
            assert!(area_tri2(&vtx2xy[i0_vtx], &vtx2xy[tris[itri0].v[inotri0]], &vtx2xy[i1_vtx]) > 1.0e-20_f64.as_());
            assert!(area_tri2(&vtx2xy[i0_vtx], &vtx2xy[i1_vtx], &vtx2xy[tris[itri0].v[inotri1]]) > 1.0e-20_f64.as_());
//            std::cout << ratio << std::endl;
            if ratio < 1.0e-20_f64.as_() {
                panic!();
            } else if ratio > 1.0_f64.as_() - 1.0e-10_f64.as_() {
                panic!();
            } else {
                let ied0 = 3 - inotri0 - inotri1;
                assert!(tris[itri0].s[ied0] < tris.len());
                /*
                # if !defined(NDEBUG)
                const unsigned
                int
                itri1 = aTri[itri0].s2[ied0];
                const unsigned
                int
                ied1 = FindAdjEdgeIndex(aTri[itri0], ied0, aTri);
                assert(aTri[itri1].s2[ied1] == itri0);
                # endif
                 */
                let res = flip_edge(itri0, ied0, vtx2tri, tris);
//        std::cout << itri0 << " " << ied0 << " " << ratio << " " << res << std::endl;
//        continue;
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
    tri_vtx: &mut Vec<DynamicTriangle>,
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
            vtx_xy[ipo1] = vtx_xy_tmp[ipo].clone();
        }
    }
    for itri in 0..tri_vtx.len() {
        for ifatri in 0..3 {
            let ipo = tri_vtx[itri].v[ifatri];
            assert_ne!(map_po_del[ipo], usize::MAX);
            tri_vtx[itri].v[ifatri] = map_po_del[ipo];
            vtx_tri[ipo].e = itri;
            vtx_tri[ipo].d = ifatri;
        }
    }
}

pub fn meshing_single_connected_shape2(
    vtx_tri: &mut Vec<DynamicVertex>,
    vtx_xy: &mut Vec<nalgebra::Vector2<f32>>,
    tri_vtx: &mut Vec<DynamicTriangle>,
    loop_vtx_idx: &Vec<usize>,
    loop_vtx: &Vec<usize>)
{
    use crate::topology::{
        assert_dynamic_triangle_mesh,
        assert_dynamic_triangles,
        find_edge_by_looking_all_triangles,
        flag_connected,
        delete_tri_flag,
    };
    let mut point_idx_to_delete = Vec::<usize>::new();
    {
        let npo = vtx_xy.len();
        point_idx_to_delete.push(npo + 0);
        point_idx_to_delete.push(npo + 1);
        point_idx_to_delete.push(npo + 2);
    }
    meshing_initialize(tri_vtx, vtx_tri, vtx_xy);
    #[cfg(debug_assertions)]
        {
            assert_dynamic_triangles(&tri_vtx);
            assert_dynamic_triangle_mesh(&vtx_tri, &tri_vtx);
        }
    for iloop in 0..loop_vtx_idx.len() - 1 {
        let nvtx = loop_vtx_idx[iloop + 1] - loop_vtx_idx[iloop];
        for iivtx in loop_vtx_idx[iloop]..loop_vtx_idx[iloop + 1] {
            let ivtx0 = loop_vtx[loop_vtx_idx[iloop] + (iivtx + 0) % nvtx];
            let ivtx1 = loop_vtx[loop_vtx_idx[iloop] + (iivtx + 1) % nvtx];
            enforce_edge(vtx_tri, tri_vtx,
                         ivtx0, ivtx1, &vtx_xy);
        }
    }
    {
        let mut aflg = vec!(0; tri_vtx.len());
        let mut itri0_ker = usize::MAX;
        let mut iedtri = 0;
        find_edge_by_looking_all_triangles(
            &mut itri0_ker, &mut iedtri,
            loop_vtx[0], loop_vtx[1], &tri_vtx);
        assert!(itri0_ker < tri_vtx.len());
        flag_connected(
            &mut aflg,
            &tri_vtx, itri0_ker, 1);
        delete_tri_flag(tri_vtx, &mut aflg, 0);
    }
    delete_unreferenced_points(
        vtx_xy, vtx_tri, tri_vtx,
        &point_idx_to_delete);
    #[cfg(debug_assertions)]
        {
            assert_dynamic_triangles(&tri_vtx);
            assert_dynamic_triangle_mesh(&vtx_tri, &tri_vtx);
        }
}

fn laplacian_mesh_smoothing_around_point<T>(
    vtx_xy: &mut Vec<nalgebra::Vector2<T>>,
    ipoin: usize,
    vtx_tri: &Vec<DynamicVertex>,
    tri_vtx: &Vec<DynamicTriangle>) -> bool
    where T: num_traits::Float + Clone + std::fmt::Debug + std::ops::AddAssign +
    std::ops::DivAssign + 'static + Copy,
          usize: num_traits::AsPrimitive<T>
{
    use crate::topology::move_ccw;
    assert_eq!(vtx_xy.len(), vtx_tri.len());
    let mut itri0 = vtx_tri[ipoin].e;
    let mut ino0 = vtx_tri[ipoin].d;
    let mut vec_delta = vtx_xy[ipoin].clone();
    let mut ntri_around: usize = 1;
    loop { // counter-clock wise
        assert!(itri0 < tri_vtx.len() && ino0 < 3 && tri_vtx[itri0].v[ino0] == ipoin);
        vec_delta += vtx_xy[tri_vtx[itri0].v[(ino0 + 1) % 3]];
        ntri_around += 1;
        if !move_ccw(&mut itri0, &mut ino0, usize::MAX, tri_vtx) { return false; }
        if itri0 == vtx_tri[ipoin].e { break; }
    }
    vtx_xy[ipoin] = vec_delta / ntri_around.as_();
    return true;
}

pub fn meshing_inside<T>(
    vtx2tri: &mut Vec<DynamicVertex>,
    tris: &mut Vec<DynamicTriangle>,
    vtx2xy: &mut Vec<nalgebra::Vector2<T>>,
    vtx2flag: &mut Vec<usize>,
    tri2flag: &mut Vec<usize>,
    num_vtx_fix: usize,
    nflgpnt_offset: usize,
    target_len: T)
    where T: num_traits::Float + std::ops::AddAssign + std::ops::DivAssign + std::ops::MulAssign<f64> +
    'static + std::fmt::Debug + std::default::Default,
          f64: AsPrimitive<T>,
          usize: AsPrimitive<T>
{
    use crate::topology::insert_a_point_inside_an_element;
    assert_eq!(vtx2xy.len(), vtx2tri.len());
    assert_eq!(vtx2flag.len(), vtx2tri.len());
    assert_eq!(tri2flag.len(), tris.len());

    let mut ratio = 3_f64.as_();
    loop {
        let mut nadd = 0;
        for itri in 0..tris.len() {
            let area = area_tri2(&vtx2xy[tris[itri].v[0]],
                                 &vtx2xy[tris[itri].v[1]],
                                 &vtx2xy[tris[itri].v[2]]);
            let _pcnt: [T; 2] = [
                (vtx2xy[tris[itri].v[0]][0] + vtx2xy[tris[itri].v[1]][0] + vtx2xy[tris[itri].v[2]][0]) / 3_f64.as_(),
                (vtx2xy[tris[itri].v[0]][1] + vtx2xy[tris[itri].v[1]][1] + vtx2xy[tris[itri].v[2]][1]) / 3_f64.as_()
            ];
            let len2 = target_len; // len * mesh_density.edgeLengthRatio(pcnt[0], pcnt[1]); //
            if area < len2 * len2 * ratio { continue; }
            let ipo0 = vtx2tri.len();
            vtx2tri.resize(vtx2tri.len() + 1, DynamicVertex { e: 0, d: 0 });
            vtx2xy.resize(vtx2xy.len() + 1, Default::default());
            vtx2xy[ipo0].x = (vtx2xy[tris[itri].v[0]].x + vtx2xy[tris[itri].v[1]].x + vtx2xy[tris[itri].v[2]].x) / 3_f64.as_();
            vtx2xy[ipo0].y = (vtx2xy[tris[itri].v[0]].y + vtx2xy[tris[itri].v[1]].y + vtx2xy[tris[itri].v[2]].y) / 3_f64.as_();
            insert_a_point_inside_an_element(ipo0, itri, vtx2tri, tris);
            let iflgtri = tri2flag[itri];
            tri2flag.push(iflgtri);
            tri2flag.push(iflgtri);
            vtx2flag.push(iflgtri + nflgpnt_offset);
            delaunay_around_point(ipo0, vtx2tri, tris, vtx2xy);
            nadd += 1;
        }
        for ip in num_vtx_fix..vtx2xy.len() {
            laplacian_mesh_smoothing_around_point(
                vtx2xy,
                ip,
                vtx2tri, tris);
        }
        if nadd != 0 { ratio *= 0.8; } else { ratio *= 0.5; }
        if ratio < 0.65.as_() {
            break;
        }
    }

    for ip in num_vtx_fix..vtx2xy.len() {
        laplacian_mesh_smoothing_around_point(
            vtx2xy,
            ip,
            vtx2tri, tris);
        delaunay_around_point(
            ip,
            vtx2tri, tris, vtx2xy);
    }
}

// --------------------------

#[test]
fn test_square() {
    use crate::topology::{
        assert_dynamic_triangle_mesh,
        assert_dynamic_triangles,
    };
    let loop2idx = vec!(0, 4);
    let idx2vtx = vec!(0, 1, 2, 3);
    type Vec2 = nalgebra::Vector2<f32>;
    let mut vtx2xy = Vec::<Vec2>::new();
    {
        vtx2xy.push(Vec2::new(-1.0, -1.0));
        vtx2xy.push(Vec2::new(1.0, -1.0));
        vtx2xy.push(Vec2::new(1.0, 1.0));
        vtx2xy.push(Vec2::new(-1.0, 1.0));
    }
    let mut tri2pnt = Vec::<DynamicTriangle>::new();
    let mut pnt2tri = Vec::<DynamicVertex>::new();
    meshing_single_connected_shape2(
        &mut pnt2tri, &mut vtx2xy, &mut tri2pnt,
        &loop2idx, &idx2vtx);
    assert_dynamic_triangles(&tri2pnt);
    assert_dynamic_triangle_mesh(&pnt2tri, &tri2pnt);
    assert_eq!(pnt2tri.len(), 4);
    assert_eq!(vtx2xy.len(), 4);
    assert_eq!(tri2pnt.len(), 2);
}

