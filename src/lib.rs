use crate::topology::DynamicTriangle;

pub mod mesher2;
pub mod topology;
mod geometry2;

pub fn array_from_2d_dynamic_triangle_mesh(
    dtris: &Vec<DynamicTriangle>,
    dvtxs: &Vec<nalgebra::Vector2<f64>>) -> (Vec<usize>, Vec<f64>) {
    let tri2vtx = {
        let mut tri2vtx = Vec::<usize>::new();
        tri2vtx.reserve(dtris.len() * 3);
        for tri in dtris.iter() {
            tri2vtx.push(tri.v[0]);
            tri2vtx.push(tri.v[1]);
            tri2vtx.push(tri.v[2]);
        }
        tri2vtx
    };
    let vtx2xy = {
        let mut vtx2xy = Vec::<f64>::new();
        vtx2xy.reserve(dvtxs.len() * 2);
        for vtx in dvtxs.iter() {
            vtx2xy.push(vtx.x);
            vtx2xy.push(vtx.y);
        }
        vtx2xy
    };
    (tri2vtx, vtx2xy)
}

pub fn area_of_triangles(
    dtris: &Vec<DynamicTriangle>,
    dvtxs: &[nalgebra::Vector2<f64>]) -> Vec<f64> {
    let mut areas = Vec::<f64>::new();
    areas.reserve(dtris.len());
    for tri in dtris {
        let p0 = &dvtxs[tri.v[0]];
        let p1 = &dvtxs[tri.v[1]];
        let p2 = &dvtxs[tri.v[2]];
        areas.push(crate::geometry2::area_tri2(p0,p1,p2));
    }
    areas
}