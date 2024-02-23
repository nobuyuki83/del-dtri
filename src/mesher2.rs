//! 2D dynamic mesh editing utility code

use num_traits::AsPrimitive;
use crate::topology::{DynamicTriangle, DynamicVertex};

#[allow(clippy::identity_op)]
pub fn meshing_single_connected_shape2(
    vtx2xy: &mut Vec<nalgebra::Vector2<f32>>,
    loop2idx: &[usize],
    idx2vtx: &[usize]) -> (Vec<DynamicTriangle>, Vec<DynamicVertex>)
{

    let mut point_idx_to_delete = Vec::<usize>::new();
    {
        let npo = vtx2xy.len();
        point_idx_to_delete.push(npo + 0);
        point_idx_to_delete.push(npo + 1);
        point_idx_to_delete.push(npo + 2);
    }
    let (mut tri2vtx, mut vtx2tri)
        = crate::trimesh2_util::meshing_initialize(vtx2xy);
    // crate::trimesh2_util::save_tri_mesh("hoge.obj", &tri_vtx, vtx_xy);
    debug_assert!(crate::topology::check_dynamic_triangle_mesh_topology(&vtx2tri, &tri2vtx));
    for iloop in 0..loop2idx.len() - 1 {
        let nvtx = loop2idx[iloop + 1] - loop2idx[iloop];
        for iivtx in loop2idx[iloop]..loop2idx[iloop + 1] {
            let ivtx0 = idx2vtx[loop2idx[iloop] + (iivtx + 0) % nvtx];
            let ivtx1 = idx2vtx[loop2idx[iloop] + (iivtx + 1) % nvtx];
            crate::trimesh2_util::enforce_edge(&mut vtx2tri, &mut tri2vtx,
                                               ivtx0, ivtx1, vtx2xy);
        }
    }
    {
        let mut aflg = vec!(0; tri2vtx.len());
        let mut itri0_ker = usize::MAX;
        let mut iedtri = 0;
        crate::topology::find_edge_by_looking_all_triangles(
            &mut itri0_ker, &mut iedtri,
            idx2vtx[0], idx2vtx[1], &tri2vtx);
        assert!(itri0_ker < tri2vtx.len());
        crate::topology::flag_connected(
            &mut aflg,
            &tri2vtx, itri0_ker, 1);
        crate::topology::delete_tri_flag(&mut tri2vtx, &mut aflg, 0);
    }
    crate::trimesh2_util::delete_unreferenced_points(
        vtx2xy, &mut vtx2tri, &mut tri2vtx,
        &point_idx_to_delete);
    debug_assert!(
        crate::topology::check_dynamic_triangle_mesh_topology(&vtx2tri, &tri2vtx));
    (tri2vtx, vtx2tri)
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
    where T: num_traits::Float + std::ops::AddAssign + std::ops::DivAssign+
    'static + std::fmt::Debug + std::default::Default + std::ops::MulAssign,
          f64: AsPrimitive<T>,
          usize: AsPrimitive<T>
{
    use crate::topology::insert_a_point_inside_an_element;
    assert_eq!(vtx2xy.len(), vtx2tri.len());
    assert_eq!(vtx2flag.len(), vtx2tri.len());
    assert_eq!(tri2flag.len(), tris.len());

    let mut ratio: T = 3_f64.as_();
    loop {
        let mut nadd = 0;
        for itri in 0..tris.len() {
            let area = del_geo::tri2::area(
                &vtx2xy[tris[itri].v[0]],
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
            crate::trimesh2_util::delaunay_around_point(ipo0, vtx2tri, tris, vtx2xy);
            nadd += 1;
        }
        for ip in num_vtx_fix..vtx2xy.len() {
            crate::trimesh2_util::laplacian_mesh_smoothing_around_point(
                vtx2xy,
                ip,
                vtx2tri, tris);
        }
        if nadd != 0 { ratio *= 0.8_f64.as_(); } else { ratio *= 0.5_f64.as_(); }
        if ratio < 0.65.as_() {
            break;
        }
    }

    for ip in num_vtx_fix..vtx2xy.len() {
        crate::trimesh2_util::laplacian_mesh_smoothing_around_point(
            vtx2xy,
            ip,
            vtx2tri, tris);
        crate::trimesh2_util::delaunay_around_point(
            ip,
            vtx2tri, tris, vtx2xy);
    }
}

// --------------------------

#[test]
fn test_square() {
    use crate::topology::check_dynamic_triangle_mesh_topology;
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
    let (tri2pnt, pnt2tri) = meshing_single_connected_shape2(
        &mut vtx2xy,
        &loop2idx, &idx2vtx);
    check_dynamic_triangle_mesh_topology(&pnt2tri, &tri2pnt);
    assert_eq!(pnt2tri.len(), 4);
    assert_eq!(vtx2xy.len(), 4);
    assert_eq!(tri2pnt.len(), 2);
}

#[test]
fn test_resampled() {
    let vtx2xy_in: Vec<f32> = vec!(
        0.0, 0.0,
        1.0, 0.0,
        1.0, 0.1,
        0.0, 0.1);
    let resolution_edge = 0.11;
    let num_vtx = vtx2xy_in.len() / 2;
    let mut loop2idx = vec!(0, num_vtx);
    let mut idx2vtx = Vec::<usize>::from_iter(0..num_vtx);
    type Vec2 = nalgebra::Vector2<f32>;
    let mut vtx2xy = vtx2xy_in.chunks(2).map(|v| Vec2::new(v[0],v[1]) ).collect();
    //
    if resolution_edge > 0. { // resample edge edge
        del_msh::polyloop::resample_multiple_loops_remain_original_vtxs(
            &mut loop2idx, &mut idx2vtx, &mut vtx2xy, resolution_edge);
    }
    //
    let (mut _tri2vtx, mut _vtx2tri)
        = crate::mesher2::meshing_single_connected_shape2(
        &mut vtx2xy,
        &loop2idx, &idx2vtx);
}