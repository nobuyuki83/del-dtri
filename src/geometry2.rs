use num_traits::AsPrimitive;

pub fn bounding_box2<VEC>(
    vtx_xy: &mut [VEC])
    -> (VEC, VEC)
    where VEC: Copy + std::ops::IndexMut<usize>,
          <VEC as std::ops::Index<usize>>::Output: PartialOrd + Sized + Copy
{
    let (mut vmin, mut vmax) = (vtx_xy[0], vtx_xy[0]);
    for xy in vtx_xy.iter().skip(1) {
        if xy[0] < vmin[0] { vmin[0] = xy[0]; }
        if xy[0] > vmax[0] { vmax[0] = xy[0]; }
        if xy[1] < vmin[1] { vmin[1] = xy[1]; }
        if xy[1] > vmax[1] { vmax[1] = xy[1]; }
    }
    (vmin, vmax)
}

pub fn area_tri2<T>(
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
    where T: num_traits::real::Real + 'static + Copy
{
    (v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[1] - v2[1]) * (v1[1] - v2[1])
}

pub fn det_delaunay<T>(
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

    let qradius = squared_distance(&out_center, p0);
    let qdistance = squared_distance(&out_center, p3);

//	assert( fabs( qradius - SquareLength(out_center,p1) ) < 1.0e-10*qradius );
//	assert( fabs( qradius - SquareLength(out_center,p2) ) < 1.0e-10*qradius );

    let tol: T = 1.0e-20.as_();
    if qdistance > qradius * (1_f64.as_() + tol) { 2 }    // outside the circumcircle
    else if qdistance < qradius * (1_f64.as_() - tol) { 0 }    // inside the circumcircle
    else { 1 }    // on the circumcircle
}