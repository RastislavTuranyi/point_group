pub mod geometry;
pub mod utils;
//pub mod point_groups;

use crate::geometry::*;

/// Finds the point group of a single molecule
/// 
/// # Example
/// 
/// ```
/// assert!(true)
/// ```
pub fn find_point_group(coordinates: Vec<Point>, precision: u16, tolerance: f64) -> String {
    if coordinates.len() == 1 {
        return "spherical".to_owned()
    } else if coordinates.len() == 2 || is_linear(&coordinates, tolerance) {
        if has_centre_of_inversion(&coordinates, tolerance) {
            return "D inf h".to_owned()
        } else {
            return "C inf v".to_owned()
        }
    } 
    
    let centre_of_symmetry = get_centre_of_symmetry(&coordinates);

    let mut rotation_axes = get_rotation_axes(&coordinates, centre_of_symmetry, precision, tolerance);
    let mut rotations = compute_rotations(&rotation_axes, &coordinates, tolerance);

    let planar = has_planar_element(&coordinates, centre_of_symmetry, tolerance);
    match planar {
        Planarity::Planar(plane) | Planarity::PlanarElement(plane) => {
            let line = Line::new(centre_of_symmetry, plane.normal);
            let rot = compute_rotations(&vec![line], &coordinates, tolerance)[0];

            rotation_axes.push(line);
            rotations.push(rot);
        },
        Planarity::None => ()
    }

    println!("PLANAR = {:?}", &planar);

    let (rotation_axes, rotations) = filter_rotations(rotation_axes, rotations, tolerance);

    let highest_order_rotation = match rotations.iter().max() {
        Some(value) => *value,
        None => 1
    };

    //println!("{:?}", &rotation_axes);
    println!("");
    println!("{}, {:?}, highest = {}", rotations.len(), rotations, highest_order_rotation);
    for (rot, axis) in rotations.iter().zip(&rotation_axes) {
        if rot >= &2 {
            println!("{} - {:?}", rot, axis);
        } else if axis.approx_eq(&Line { point: Point { x: 7.1446236518, y: 9.656128203933331, z: 4.7408198414000005 }, vector: Vector { x: -0.3112557392670835, y: 0.5601211141467458, z: 1.8173851621853414 } }, tolerance) {
            println!("SIMILAR {:?}", axis)
        };

        /* println!("{:?} {:?}", deviation_from_rotation_symmetry(&coordinates, *axis, 2.0).iter().fold(0.0, |max, val| if val > &max {*val} else {max}), axis);

        if axis.vector.is_approx_k_multiple(Vector { x: 0.8808450198738458, y: 0.4532289546188045, z: -0.13673172879251713 }, tolerance) {
            print!("THIS ONE!!!!!!!!!!!!!!!!!!!!!!!");
        } */

        /* for r in rotation_axes[1..].iter() {
            println!("{}", axis.approx_eq(r, tolerance))
        } */
    }
    
    // Check if the compound has a high symmetry, and if it has no rotational symmetry
    if highest_order_rotation >= 3 {
        let number_of_highest_rotations = rotations.iter()
            .filter(|x| **x == highest_order_rotation)
            .count();
        
        if number_of_highest_rotations > 1 {
            if !has_centre_of_inversion(&coordinates, tolerance) {
                return "T d".to_owned()
            }

            let c5 = rotations.iter().find(|x| **x == 5);

            match c5 {
                Some(_) => return "I h".to_owned(),
                None => return "O h".to_owned()
            }
        }
    } else if highest_order_rotation == 1 {
        match find_reflection_plane(&coordinates, &centre_of_symmetry, precision, tolerance) {
            Some(_plane) => return "C s".to_owned(),
            None => ()
        }

        if has_centre_of_inversion(&coordinates, tolerance) {
            return "C i".to_owned()
        }

        return "C 1".to_owned()
    }

    let number_of_c2 = rotations.iter().filter(|x| x == &&2_u32).count() as u32;

    if number_of_c2 == highest_order_rotation || (highest_order_rotation == 2 && number_of_c2 == highest_order_rotation+1) {
        let principal_axis = rotation_axes[rotations.iter()
            .position(|x| x == &highest_order_rotation)
            .expect("We are searching for a value we have retrieved from here earlier")];

        let sigma_h = get_horizontal_plane(&coordinates, &principal_axis, &centre_of_symmetry, tolerance);

        match sigma_h {
            Some(_i) => return format!("D {} h", highest_order_rotation),
            None => ()
        };

        let sigma_v = if highest_order_rotation == 2 {
            get_d2d_vertical_planes(&coordinates, &rotation_axes, precision, planar, tolerance)
        } else {
            get_vertical_planes(&coordinates, &principal_axis, precision, planar, tolerance)
        };
        println!("{:?} {}", &sigma_v, sigma_v.len());

        if (sigma_v.len() as u32) == highest_order_rotation {
            return format!("D {} d", highest_order_rotation)
        } else {
            return format!("D {}", highest_order_rotation)
        }
    } else {
        let principal_axis = rotation_axes[rotations.iter()
            .position(|x| x == &highest_order_rotation)
            .expect("We are searching for a value we have retrieved from here earlier")];
        println!("{:?}", &principal_axis);

        let sigma_h = get_horizontal_plane(&coordinates, &principal_axis, &centre_of_symmetry, tolerance);
        println!("sigma h {:?}", &sigma_h);

        match sigma_h {
            Some(_i) => return format!("C {} h", highest_order_rotation),
            None => ()
        };

        let sigma_v = get_vertical_planes(&coordinates, &principal_axis, precision, planar, tolerance);
        println!("sigma v {:?}", &sigma_v);

        if (sigma_v.len() as u32) == highest_order_rotation {
            return format!("C {} v", highest_order_rotation)
        };

        let plane = Plane{point: centre_of_symmetry.clone(), normal: principal_axis.vector.clone()};
        
        if is_s2n_present(&coordinates, &rotation_axes, highest_order_rotation.clone() as f64, plane, tolerance) {
            return format!("S {}", 2 * highest_order_rotation)
        } else {
            return format!("C {}", highest_order_rotation)
        }
    }
}



fn is_linear(coordinates: &Vec<Point>, tolerance: f64) -> bool {
    let line = Line::new(coordinates[0].clone(), Vector::from_two_points(coordinates[0], coordinates[1]));

    // Determine wheter each subsequent point lies on the line, i.e. if vectors are k multiple of one another
    for coords in &coordinates[2..] {
        if !line.approx_has_point(*coords, tolerance) {
            return false
        }
    }

    return true
}


fn has_centre_of_inversion(coordinates: &Vec<Point>, tolerance: f64) -> bool {
    let centre_of_symmetry = get_centre_of_symmetry(coordinates);

    for point in coordinates {
        let inverted = point.invert_around(centre_of_symmetry);

        if !is_equivalent_point(&coordinates, inverted, tolerance) {
            return false
        }
    }

    return true
}


fn get_rotation_axes(coordinates: &Vec<Point>, centre_of_symmetry: Point, precision: u16, tolerance: f64) -> Vec<Line> {
    let length = coordinates.len();

    // Create a vector to hold the axes with at least the exactly required size (hopefully)
    let mut axes = if precision == 1 {
        Vec::<Line>::with_capacity(length)
    } else if precision == 2 {
        Vec::<Line>::with_capacity(combinations(length, 2) + length)
    } else if precision == 3 {
        Vec::<Line>::with_capacity(length + combinations(length, 2) + combinations(length, 3))
    } else if precision == 4 {
        Vec::<Line>::with_capacity(length + combinations(length, 2) + 10)
    } else if precision == 5 {
        Vec::<Line>::with_capacity(length + combinations(length, 2) + combinations(length, 3) + 10)
    } else {
        Vec::with_capacity(0)
    };

    let coords = if precision > 1 {
        create_pseudo_points(coordinates, precision, tolerance)
    } else {
        coordinates.clone()
    };

    // Create lines between the centre of symmetry and each (pseudo) atom
    for point in coords {
        let line = match Line::from_points_tolerance_normalised(centre_of_symmetry, point, tolerance) {
            Some(line) => line.normalise(),
            None => continue
        };

        axes.push(line)
    }
    
    return axes
}



fn _is_equivalent_axis(axes: &Vec<Line>, line: Line, tolerance: f64) -> bool {
    for ax in axes {
        if ax.approx_eq(&line, tolerance) {
            return true
        }
    }

    return false
}



fn _get_rotation_axes_thorough(coordinates: &Vec<Point>, precision: u16, tolerance: f64) -> Vec<Line> {
    let length = coordinates.len();

    // Create a vector to hold the axes with at least the exactly required size (hopefully)
    let mut axes = if precision == 1 {
        Vec::<Line>::with_capacity(combinations(length, 2))
    } else if precision == 2 {
        Vec::<Line>::with_capacity(combinations(length + combinations(length, 2), 2))
    } else if precision == 3 {
        Vec::<Line>::with_capacity(combinations(length + combinations(length, 2) + combinations(length, 3), 2))
    } else if precision == 4 {
        Vec::<Line>::with_capacity(combinations(length + combinations(length, 2), 2) + 10)
    } else if precision == 5 {
        Vec::<Line>::with_capacity(combinations(length + combinations(length, 2) + combinations(length, 3), 2) + 10)
    } else {
        Vec::with_capacity(0)
    };

    let coords = if precision > 1 {
        create_pseudo_points(coordinates, precision, tolerance)
    } else {
        coordinates.clone()
    };

    let length = coords.len();
    // Create lines between all pairs of atoms
    let end = length - 1;
    for (i, outer_coords) in coords.iter().enumerate() {
        if i == end {
            break
        }

        'inner_coords: for inner_coords in coords[(i+1)..].iter() {

            let line = match Line::from_points(*outer_coords, *inner_coords) {
                Some(line) => line,
                None => continue 'inner_coords
            };

            // Do not add this axis if an identical one already exists
            for ax in &axes {
                if ax.approx_eq(&line, tolerance) {
                    continue 'inner_coords;
                }
            }

            axes.push(line)
        }
    }

    axes
}



fn create_pseudo_points(coordinates: &Vec<Point>, precision: u16, tolerance: f64) -> Vec<Point> {
    let length = coordinates.len();

    let mut points = if precision == 2 {
        Vec::with_capacity(length + combinations(length, 2))
    } else if precision == 3 {
        Vec::with_capacity(length + combinations(length, 2) + combinations(length, 3))
    } else if precision == 4 {
        Vec::with_capacity(length + combinations(length, 2) + 10)
    } else if precision == 5 {
        Vec::with_capacity(length + combinations(length, 2) + combinations(length, 3) + 10)
    } else {
        Vec::with_capacity(0)
    };

    
    if precision == 2 {
        for (i, outer_coords) in coordinates.iter().enumerate() {
            points.push(*outer_coords);
            
            for inner_coords in coordinates[(i+1)..].iter() {
                points.push((*outer_coords + *inner_coords) / 2.0)
            }
        }
    } else if precision == 3 {
        for (i, outer) in coordinates.iter().enumerate() {
            points.push(*outer);
            
            for (j, middle) in coordinates[(i+1)..].iter().enumerate() {
                points.push((*outer + *middle) / 2.0);
            
                for inner in coordinates[(i+j+2)..].iter() {
                    points.push((*outer + *middle + *inner) / 3.0);
                }
            }
        }
    } else if precision == 4 || precision == 5 {
        let mut surfaces = Vec::with_capacity(10);

        for (i, outer) in coordinates.iter().enumerate() {
            points.push(*outer);
            
            for (j, middle) in coordinates[(i+1)..].iter().enumerate() {
                points.push((*outer + *middle) / 2.0);
            
                'inner: for (k, inner) in coordinates[(i+j+2)..].iter().enumerate() {
                     if precision == 5 {
                        points.push((*outer + *middle + *inner) / 3.0)
                     }

                    let plane = match Plane::from_three_points(*outer, *middle, *inner) {
                        Some(plane) => plane,
                        None => continue 'inner
                    };

                    if !surfaces.contains(&plane) {
                        let mut new = Vec::with_capacity(8);
                        for point in coordinates[(i+j+k+3)..].iter() {
                            if point.is_approx_on_plane(plane, tolerance) {
                                new.push(point)
                            }
                        }

                        if !new.is_empty() {
                            let mut sum = *outer + *inner + *middle;
                            for point in &new {
                                sum += **point;
                            }
                            points.push(sum / (new.len() as f64 + 3.0));
                            surfaces.push(plane)
                        }
                    }
                    
                }
            }
        }
    }

    return points
}


pub fn is_planar(coordinates: &Vec<Point>, centre_of_symmetry: Point, tolerance: f64) -> Option<Plane> {
    let plane = match get_valid_plane(&coordinates, centre_of_symmetry, tolerance) {
        Some(plane) => plane,
        None => return None
    };

    for point in coordinates[3..].iter() {
        if !plane.approx_has_point(*point, tolerance) {
            return None
        }
    }

    return Some(plane)
}



pub fn has_planar_element(coordinates: &Vec<Point>, centre_of_symmetry: Point, tolerance: f64) -> Planarity {
    let plane = match get_valid_plane(&coordinates, centre_of_symmetry, tolerance) {
        Some(plane) => plane,
        None => return Planarity::None
    };

    let (plane, _, length) = find_best_plane(&coordinates, plane, tolerance);
    
    if length == coordinates.len() {
        return Planarity::Planar(plane)
    } else if length > 4 {
        return Planarity::PlanarElement(plane)
    } else {
        return Planarity::None
    }
}


#[derive(Debug, Copy, Clone)]
pub enum Planarity {
    Planar(Plane),
    PlanarElement(Plane),
    None
}



fn get_valid_plane(points: &Vec<Point>, central_point: Point, tolerance: f64) -> Option<Plane> {
    for (i, point2) in points.iter().enumerate() {
        for point3 in points[i+1..].iter() {
            match Plane::from_three_points_tolerance(central_point, *point2, *point3, tolerance) {
                Some(plane) => return Some(plane),
                None => continue
            }
        }
    }

    return None
}



fn find_points_on_plane(points: &Vec<Point>, plane: Plane, tolerance: f64) -> Vec<&Point> {
    points.iter().filter(|point| plane.approx_has_point(**point, tolerance)).collect::<Vec<_>>()
}


fn find_best_plane(coordinates: &Vec<Point>, plane: Plane, tolerance: f64) -> (Plane, Vec<&Point>, usize) {
    let mut points_on_plane = find_points_on_plane(&coordinates, plane, tolerance);
    let mut length = points_on_plane.len();
    let mut new_length;

    loop {
        let plane = match Plane::from_multiple_points(&points_on_plane, tolerance) {
            Some(plane) => plane,
            None => return (plane, points_on_plane, length)
        };

        points_on_plane = find_points_on_plane(&coordinates, plane, tolerance);
        new_length = points_on_plane.len();

        if new_length <= length {
            return (plane, points_on_plane, new_length)
        }

        length = new_length;
    }
}


pub fn compute_rotations(axes: &Vec<Line>, coordinates: &Vec<Point>, tolerance: f64) -> Vec<u32> {
    let mut rotations = Vec::<u32>::with_capacity(axes.len());

    for axis in axes {
        let mut rot = 0.0;

        for point in coordinates {
            if rot == 0.0 {
                if point.distance(*axis) < tolerance {
                    continue
                }

                rot = determine_baseline_rotation(&coordinates, point, axis, tolerance);

                if rot == 1.0 {
                    break
                }
            } else {
                if !is_rotation_viable(&coordinates, point, axis, rot, tolerance) {
                    rot = recompute_baseline_rotation(coordinates, point, axis, rot, tolerance);

                    if rot == 1.0 {
                    break
                    }
                }
            }
        }
            
        rotations.push(rot as u32);
    }

    return rotations
}


fn filter_rotations(axes: Vec<Line>, rotations: Vec<u32>, tolerance: f64) -> (Vec<Line>, Vec<u32>) {
    let length = rotations.len();
    let mut nonzero_axes = Vec::with_capacity(length);
    let mut nonzero_rots = Vec::with_capacity(length);

    'outer: for (i, (axis, rot)) in axes.iter().zip(rotations.iter()).enumerate() {
        if *rot == 1 {
            continue;
        }

        for (comparison_axis, comparison_rot) in axes[i+1..].iter().zip(rotations[i+1..].iter()) {
            if *rot == 1 {
                continue;
            }

            if axis.fast_approx_eq(*comparison_axis, tolerance) && rot <= comparison_rot {
                continue 'outer
            }
        }

        nonzero_axes.push(*axis);
        nonzero_rots.push(*rot);
    }

    return (nonzero_axes, nonzero_rots)
}


fn determine_baseline_rotation(coordinates: &Vec<Point>, point: &Point, axis: &Line, tolerance: f64) -> f64 {
    // Finds the rotation with the smallest rotation (Cn) that preserves the symmetry of the given object
    // This function is meant to be used only on the first atom of a compound since it takes into account all n
    // from 1 to 8.
    if is_rotation_viable(&coordinates, &point, &axis, 2.0, tolerance) {
        if is_rotation_viable(&coordinates, &point, &axis, 4.0, tolerance) {
            if is_rotation_viable(&coordinates, &point, &axis, 8.0, tolerance) {
                return 8.0
            } else {
                return 4.0
            }
        } else if is_rotation_viable(&coordinates, &point, &axis, 6.0, tolerance) {
            return 6.0
        } else {
            return 2.0
        }
    } else if is_rotation_viable(&coordinates, &point, &axis, 3.0, tolerance) {
        3.0
    } else if is_rotation_viable(&coordinates, &point, &axis, 5.0, tolerance) {
        5.0
    } else if is_rotation_viable(&coordinates, &point, &axis, 7.0, tolerance) {
        7.0
    } else {
        1.0
    }
}


fn recompute_baseline_rotation(coordinates: &Vec<Point>, point: &Point, axis: &Line, rot: f64, tolerance: f64) -> f64 {
    // Computes the smallest angle rotation Cn that preserves the symmetry of a given object but that is also
    // compatible with the previous smallest rotation (i.e. the smallest Cn that is a superset of the previous Cn)
    if rot == 8.0 {
        if is_rotation_viable(&coordinates, &point, &axis, 4.0, tolerance) {
            return 4.0
        } else if is_rotation_viable(&coordinates, &point, &axis, 2.0, tolerance) {
            return 2.0
        } else {
            return 1.0
        }
    } else if rot == 4.0 {
        if is_rotation_viable(&coordinates, &point, &axis, 2.0, tolerance) {
            return 2.0
        } else {
            return 1.0
        }
    } else if rot == 6.0 {
        if is_rotation_viable(&coordinates, &point, &axis, 2.0, tolerance) {
            return 2.0
        } else if is_rotation_viable(&coordinates, &point, &axis, 3.0, tolerance) {
            return 3.0
        } else {
            return 1.0
        }
    } else {
        return 1.0
    }
}


pub fn is_rotation_viable(coordinates: &Vec<Point>, point: &Point, axis: &Line, n: f64, tolerance: f64) -> bool {
    let rotated = axis.rotate_point(*point, n);
    // TODO: Check that the rotation is viable in both directions (e.g. that C3 == 3C3)
    return is_equivalent_point(&coordinates, rotated, tolerance)
}


fn is_equivalent_point(coordinates: &Vec<Point>, point: Point, tolerance: f64) -> bool {
    // Determines whether a given point is equivalent (at the same position) to any of the provided coordinates
    for coord in coordinates {
        if coord.approx_eq(&point, tolerance) {
            return true
        }
    }

    return false
}


fn combinations(n: usize, r: usize) -> usize {
    (1..=r.min(n - r)).fold(1, |acc, val| acc * (n - val + 1) / val)
}


pub fn get_centre_of_symmetry(coordinates: &Vec<Point>) -> Point {
    // Finds the centre of symmetry of a compound by averaging all the coordinates
    let mut sum = Point::new(0.0, 0.0, 0.0);

    for point in coordinates {
        sum += *point;
    }

    let len = coordinates.len() as f64;

    return sum / len
}


fn get_horizontal_plane(coordinates: &Vec<Point>, principal_axis: &Line, centre_of_symmetry: &Point, tolerance: f64) -> Option<Plane> {
    // Finds the horizontal mirror plane of the compound if one exists
    let plane = match Plane::new(centre_of_symmetry.clone(), principal_axis.vector.clone()) {
        Some(plane) => plane,
        None => return None
    };

    for point in coordinates {
        if point.is_approx_on_plane(plane, tolerance) {
            continue
        }
        let reflected = plane.reflect_point(*point);

        if !is_equivalent_point(&coordinates, reflected, tolerance) {
            return None
        }
    }

    return Some(plane)
}


fn get_vertical_planes(coordinates: &Vec<Point>, principal_axis: &Line, precision: u16, 
        planar: Planarity, tolerance: f64) -> Vec<Plane> {
    let length = coordinates.len();
    
    let mut planes = if precision == 1 {
        Vec::<Plane>::with_capacity(combinations(length, 2))
    } else if precision >= 2 {
        Vec::<Plane>::with_capacity(combinations(length + combinations(length, 2), 2))
    } else {
        Vec::with_capacity(0)
    };

    if precision == 1 {
        'outer: for point in coordinates {
            for plane in &planes {
                if plane.has_point(*point) {
                    continue 'outer
                }
            }
            
            let vector = Vector::from_two_points(principal_axis.point, *point);
            let plane = Plane{point: *point, normal: vector.cross_product(principal_axis.vector)};

            if has_plane_reflection_symmetry(&coordinates, &plane, tolerance) {
                planes.push(plane)
            }
        }
    } else {
        let pseudo_points = create_pseudo_points(&coordinates, 2, tolerance);

        'outer: for point in pseudo_points {
            let vector = Vector::from_two_points(principal_axis.point, point);
            let plane = match Plane::new(point, vector.cross_product(principal_axis.vector)) {
                Some(plane) => plane,
                None => continue 'outer
            };

            if has_plane_reflection_symmetry(&coordinates, &plane, tolerance) {
                for existing in &planes {
                    if existing.approx_eq(&plane, tolerance)  {
                        continue 'outer
                    }
                }
                planes.push(plane)
            }
        }
    }

    match planar {
        Planarity::Planar(plane) => { 
            let plane = match Plane::new(principal_axis.point, plane.normal.cross_product(principal_axis.vector)) {
                Some(plane) => plane,
                None => return planes
            };
            for existing in &planes {
                if existing.approx_eq(&plane, tolerance) {
                    return planes
                }
            }
            planes.push(plane)
        },
        Planarity::PlanarElement(_) | Planarity::None => ()
    }

    return planes
}


fn get_d2d_vertical_planes(coordinates: &Vec<Point>, axes: &Vec<Line>, precision: u16, 
    planar: Planarity, tolerance: f64) -> Vec<Plane> {
    for axis in axes {
        let planes = get_vertical_planes(&coordinates, axis, precision, planar, tolerance);
        if planes.len() == 2 {
            return planes
        }
    }

    return Vec::with_capacity(0)
}


pub fn has_plane_reflection_symmetry(coordinates: &Vec<Point>, plane: &Plane, tolerance: f64) -> bool {
    for point in coordinates {
        if plane.approx_has_point(*point, tolerance) {
            continue
        }

        let reflected = plane.reflect_point(*point);

        if !is_equivalent_point(&coordinates, reflected, tolerance) {
            return false
        }
    }

    return true
}


fn is_s2n_present(coordinates: &Vec<Point>, axes: &Vec<Line>, highest_order_rotation: f64, plane: Plane, tolerance: f64) -> bool {
    for axis in axes {
        if is_improper_axis_viable(&coordinates, *axis, 2.0 * highest_order_rotation, plane, tolerance) {
            return true
        }
    }

    return false
}


pub fn is_improper_axis_viable(coordinates: &Vec<Point>, axis: Line, n: f64, plane: Plane, tolerance: f64) -> bool {
    for point in coordinates {
        let rotated = point.improper_rotate_around_axis(axis, n, plane);

        if !is_equivalent_point(coordinates, rotated, tolerance) {
            return false
        }
    }

    return true
}



pub fn find_reflection_plane(coordinates: &Vec<Point>, centre_of_symmetry: &Point, precision: u16, tolerance: f64) -> Option<Plane> {
    let pseudo_points = create_pseudo_points(&coordinates, precision, tolerance);

    for (i, outer) in pseudo_points.iter().enumerate() {
        for inner in pseudo_points[i..].iter() {
            let plane = match Plane::from_three_points(*centre_of_symmetry, *outer, *inner) {
                Some(plane) => plane,
                None => continue
            };

            if has_plane_reflection_symmetry(&coordinates, &plane, tolerance) {
                return Some(plane)
            }
        }
    }

    return None
}





//////////////////////////////////////////// TESTS ////////////////////////////////////////////
#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]

    use super::*;
    use rstest::*;
    use std::f64::consts::PI;

    fn octagon_coordinates() -> Vec<Point> {
        vec!(Point::new(1.0, 0.0, 0.0), Point::new(2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 0.0), Point::new(0.0, 1.0, 0.0), 
            Point::new(-2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new(-2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 0.0), 
            Point::new(0.0, -1.0, 0.0), Point::new(2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 0.0))
    }

    fn polygon_coordinates(n_vertices: f64, radius: f64, displacement: f64, z: f64) -> Vec<Point> {
        let mut coordinates = Vec::with_capacity(n_vertices as usize);

        for i in 0..n_vertices as u64 {
            let i = i as f64;
            coordinates.push(Point::new(radius * ((2.0 * i + displacement) * PI / n_vertices).sin(), radius * ((2.0 * i + displacement) * PI / n_vertices).cos(), z))
        }

        return coordinates
    }

    fn parallel_polygons(n_vertices: f64, radius: f64, displacement: f64, centre: bool) -> Vec<Point> {
        let mut c = polygon_coordinates(n_vertices, radius, displacement, 1.0);
        let mut cc = polygon_coordinates(n_vertices, radius, displacement, -1.0);
        c.append(&mut cc);

        if centre {
            c.push(Point::new(0.0, 0.0, 0.0))
        }

        return c
    }

    fn generate_Cnh_coordinates(n_vertices: f64, double: bool) -> Vec<Point> {
        let mut coordinates = polygon_coordinates(n_vertices, 1.0, 0.0, 0.0);
        let axis = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));

        if double {
            coordinates.extend(polygon_coordinates(n_vertices, 2.0, 0.0, 0.0));
            coordinates.extend(axis.rotate_coordinates(polygon_coordinates(n_vertices, 3.0, 0.0, 0.0), 
                n_vertices.powi(2), n_vertices as usize))
        } else {
            coordinates.extend(axis.rotate_coordinates(polygon_coordinates(n_vertices, 1.5, 0.0, 0.0), 
                n_vertices.powi(2), n_vertices as usize));
            coordinates.push(Point::new(0.0, 0.0, 0.0))
        }

        return coordinates
    }

    fn generate_Cnv_coordinates(n_vertices: f64) -> Vec<Point> {
        let mut coordinates = polygon_coordinates(n_vertices, 1.0, 0.0, 1.0);
        coordinates.push(Point::new(0.0, 0.0, 0.0));
        coordinates.push(Point::new(0.0, 0.0, -1.0));

        coordinates
    }

    #[rstest]
    #[case::spherical(vec!(Point::new(0.0, 0.0, 0.0)), 0, "spherical")]
    #[case::D_inf_h_3(vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)), 0, "D inf h")]
    #[case::D_inf_h_4(vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(3.0, 0.0, 0.0)), 0, "D inf h")]
    #[case::C_inf_v(vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(3.0, 0.0, 0.0)), 0, "C inf v")]
    #[case::Td(vec!(Point::new(1.0, 1.0, 1.0), Point::new(-1.0, 1.0, -1.0), Point::new(-1.0, -1.0, 1.0), Point::new(1.0, -1.0, -1.0)), 3, "T d")]
    #[case::Oh(vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0), 
        Point::new(0.0, -1.0, 0.0), Point::new(0.0, 0.0, 1.0), Point::new(0.0, 0.0, -1.0)), 3, "O h")]
    #[case::Ih(vec!(Point::new(1.0, 0.0, 0.0), 
        Point::new(1.0/5_f64.sqrt(), 2.0/5_f64.sqrt(), 0.0),
        Point::new(1.0/5_f64.sqrt(), (5.0-5_f64.sqrt())/10.0, -((5.0+5_f64.sqrt())/10.0).sqrt()),
        Point::new(1.0/5_f64.sqrt(), (-5.0-5_f64.sqrt())/10.0, ((5.0-5_f64.sqrt())/10.0).sqrt()),
        Point::new(1.0/5_f64.sqrt(), (-5.0-5_f64.sqrt())/10.0, -((5.0-5_f64.sqrt())/10.0).sqrt()),
        Point::new(1.0/5_f64.sqrt(), (5.0-5_f64.sqrt())/10.0, ((5.0+5_f64.sqrt())/10.0).sqrt()), 
        Point::new(-1.0, 0.0, 0.0),
        Point::new(-1.0/5_f64.sqrt(), -2.0/5_f64.sqrt(), 0.0),
        Point::new(-1.0/5_f64.sqrt(), (-5.0+5_f64.sqrt())/10.0, -((5.0+5_f64.sqrt())/10.0).sqrt()),
        Point::new(-1.0/5_f64.sqrt(), (5.0+5_f64.sqrt())/10.0, -((5.0-5_f64.sqrt())/10.0).sqrt()),
        Point::new(-1.0/5_f64.sqrt(), (5.0+5_f64.sqrt())/10.0, ((5.0-5_f64.sqrt())/10.0).sqrt()),
        Point::new(-1.0/5_f64.sqrt(), (-5.0+5_f64.sqrt())/10.0, ((5.0+5_f64.sqrt())/10.0).sqrt())), 2, "I h")]
    #[case::D2h(vec!(Point::new(1.0, 0.0, 0.0), Point::new(2.0, 1.0, 0.0), Point::new(2.0, -1.0, 0.0), 
        Point::new(-1.0, 0.0, 0.0), Point::new(-2.0, 1.0, 0.0), Point::new(-2.0, -1.0, 0.0)), 2, "D 2 h")]
    #[case::D3h(vec!(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 3_f64.sqrt()/3.0, 0.0), Point::new(-0.5, -3_f64.sqrt()/6.0, 0.0), 
        Point::new(0.5, -3_f64.sqrt()/6.0, 0.0)), 3, "D 3 h")]
    #[case::D4h(vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), 
        Point::new(0.0, 1.0, 0.0), Point::new(0.0, -1.0, 0.0)), 2, "D 4 h")]
    #[case::D5h(polygon_coordinates(5.0, 1.0, 0.0, 0.0), 2, "D 5 h")]
    #[case::D6h(polygon_coordinates(6.0, 1.0, 0.0, 0.0), 2, "D 6 h")]
    #[case::D6h_double(parallel_polygons(6.0, 1.0, 0.0, true), 2, "D 6 h")]
    #[case::D7h(polygon_coordinates(7.0, 1.0, 0.0, 0.0), 3, "D 7 h")]
    #[case::D8h(octagon_coordinates(), 2, "D 8 h")]
    #[case::D2d(vec!(Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new(2.0, 1.0, 0.0), Point::new(2.0, -1.0, 0.0), Point::new(-2.0, 0.0, 1.0),
        Point::new(-2.0, 0.0, -1.0)), 3, "D 2 d")]
    #[case::D3d(vec!(Point::new(0.0, 0.0, 1.0), Point::new(0.0, 3_f64.sqrt()/3.0, 1.0), Point::new(-0.5, -3_f64.sqrt()/6.0, 1.0), 
        Point::new(0.5, -3_f64.sqrt()/6.0, 1.0), Point::new(0.0, 0.0, -1.0), Point::new(0.0, -3_f64.sqrt()/3.0, -1.0), Point::new(-0.5, 3_f64.sqrt()/6.0, -1.0), 
        Point::new(0.5, 3_f64.sqrt()/6.0, -1.0)), 3, "D 3 d")]
    #[case::D4d(vec!(Point::new(1.0, 0.0, 1.0), Point::new(-1.0, 0.0, 1.0), Point::new(0.0, 1.0, 1.0), Point::new(0.0, -1.0, 1.0), 
        Point::new(2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, -1.0), Point::new(-2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, -1.0), 
        Point::new(-2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, -1.0), Point::new(2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, -1.0)), 2, "D 4 d")]
    #[case::D5d(vec!(Point::new(0.0, 1.0, 1.0), Point::new((2.0 * PI / 5.0).sin(), (2.0 * PI / 5.0).cos(), 1.0), 
        Point::new((4.0 * PI / 5.0).sin(), -(PI / 5.0).cos(), 1.0),  Point::new(-(4.0 * PI / 5.0).sin(), -(PI / 5.0).cos(), 1.0), 
        Point::new(-(2.0 * PI / 5.0).sin(), (2.0 * PI / 5.0).cos(), 1.0), Point::new(0.0, 0.0, 0.0), Point::new(0.0, -1.0, -1.0), 
        Point::new((2.0 * PI / 5.0).sin(), -(2.0 * PI / 5.0).cos(), -1.0), Point::new((4.0 * PI / 5.0).sin(), (PI / 5.0).cos(), -1.0),  
        Point::new(-(4.0 * PI / 5.0).sin(), (PI / 5.0).cos(), -1.0), Point::new(-(2.0 * PI / 5.0).sin(), -(2.0 * PI / 5.0).cos(), -1.0)), 4, "D 5 d")]
    #[case::D6d(vec!(Point::new(1.0, 0.0, 1.0), Point::new(-1.0, 0.0, 1.0), Point::new(0.5, 3_f64.sqrt()/2.0, 1.0), Point::new(0.0, 0.0, 0.0),
        Point::new(-0.5, 3_f64.sqrt()/2.0, 1.0), Point::new(-0.5, -3_f64.sqrt()/2.0, 1.0), Point::new(0.5, -3_f64.sqrt()/2.0, 1.0), 
        Point::new((PI / 6.0).cos(), (PI/ 6.0).sin(), -1.0), Point::new((PI * 3.0 / 6.0).cos(), (PI * 3.0 / 6.0).sin(), -1.0), 
        Point::new((PI * 5.0 / 6.0).cos(), (PI * 5.0 / 6.0).sin(), -1.0), Point::new((PI * 7.0 / 6.0).cos(), (PI * 7.0 / 6.0).sin(), -1.0), 
        Point::new((PI * 9.0 / 6.0).cos(), (PI * 9.0 / 6.0).sin(), -1.0), Point::new((PI * 11.0 / 6.0).cos(), (PI * 11.0 / 6.0).sin(), -1.0)), 2, "D 6 d")]
    #[case::D7d(vec!(Point::new(0.0, 1.0, 1.0), Point::new((2.0 * PI / 7.0).sin(), (2.0 * PI / 7.0).cos(), 1.0), Point::new((4.0 * PI / 7.0).sin(), (4.0 * PI / 7.0).cos(), 1.0), 
        Point::new((6.0 * PI / 7.0).sin(), (6.0 * PI / 7.0).cos(), 1.0), Point::new((8.0 * PI / 7.0).sin(), (8.0 * PI / 7.0).cos(), 1.0), 
        Point::new((10.0 * PI / 7.0).sin(), (10.0 * PI / 7.0).cos(), 1.0), Point::new((12.0 * PI / 7.0).sin(), (12.0 * PI / 7.0).cos(), 1.0),
        Point::new(0.0, -1.0, -1.0), Point::new((2.0 * PI / 7.0).sin(), -(2.0 * PI / 7.0).cos(), -1.0), Point::new((4.0 * PI / 7.0).sin(), -(4.0 * PI / 7.0).cos(), -1.0), 
        Point::new((6.0 * PI / 7.0).sin(), -(6.0 * PI / 7.0).cos(), -1.0), Point::new((8.0 * PI / 7.0).sin(), -(8.0 * PI / 7.0).cos(), -1.0), 
        Point::new((10.0 * PI / 7.0).sin(), -(10.0 * PI / 7.0).cos(), -1.0), Point::new((12.0 * PI / 7.0).sin(), -(12.0 * PI / 7.0).cos(), -1.0),
        Point::new(0.0, 0.0, 0.0)), 4, "D 7 d")]
    #[case::D8d(vec!(Point::new(1.0, 0.0, 1.0), Point::new(2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 1.0), Point::new(0.0, 1.0, 1.0), 
        Point::new(-2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 1.0), Point::new(-1.0, 0.0, 1.0), Point::new(-2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 1.0), 
        Point::new(0.0, -1.0, 1.0), Point::new(2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 1.0), Point::new(0.0, 0.0, 0.0),
        Point::new((PI / 8.0).sin(), (PI / 8.0).cos(), -1.0), Point::new((PI * 3.0 / 8.0).sin(), (PI * 3.0 / 8.0).cos(), -1.0), 
        Point::new((PI * 5.0 / 8.0).sin(), (PI * 5.0 / 8.0).cos(), -1.0), Point::new((PI * 7.0 / 8.0).sin(), (PI * 7.0 / 8.0).cos(), -1.0), 
        Point::new((PI * 9.0 / 8.0).sin(), (PI * 9.0 / 8.0).cos(), -1.0), Point::new((PI * 11.0 / 8.0).sin(), (PI * 11.0 / 8.0).cos(), -1.0), 
        Point::new((PI * 13.0 / 8.0).sin(), (PI * 13.0 / 8.0).cos(), -1.0), Point::new((PI * 15.0 / 8.0).sin(), (PI * 15.0 / 8.0).cos(), -1.0)), 2, "D 8 d")]
    /* #[case::D2_fail(vec!(Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new(0.5, 3_f64.sqrt()/2.0, 0.0), Point::new(-0.5, 3_f64.sqrt()/2.0, 0.0), 
        Point::new(-0.5, -3_f64.sqrt()/2.0, 0.0), Point::new(0.5, -3_f64.sqrt()/2.0, 0.0), Point::new(4.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), 
        Point::new(3.5, -0.75, 0.43301270189221946), Point::new(2.5, -0.75, 0.43301270189221946), Point::new(2.5, 0.75, -0.43301270189221946), 
        Point::new(3.5, 0.75, -0.43301270189221946)), 3, "D 2")] */
    #[case::D2(vec!(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 2.0, 0.0), Point::new(0.0, -2.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0),
        Point::new(0.0, 0.0, 1.0), Point::new(0.0, 0.0, -1.0), Point::new(0.5, 0.2, 2.0), Point::new(-0.5, -0.2, 2.0), Point::new(0.5, -0.2, -2.0), 
        Point::new(-0.5, 0.2, -2.0)), 2, "D 2")]
    #[case::D3(vec!(Point::new(0.0, 0.0, 1.0), Point::new(0.0, 3_f64.sqrt()/3.0, 1.0), Point::new(-0.5, -3_f64.sqrt()/6.0, 1.0), Point::new(0.5, -3_f64.sqrt()/6.0, 1.0), 
        Point::new(0.0, 0.0, -1.0), Point::new(0.40824829046386296, -0.40824829046386296, -1.0), Point::new(-0.5576775358252052, -0.1494292453613423, -1.0), 
        Point::new(0.14942924536134225, 0.5576775358252053, -1.0)), 3, "D 3")]
    #[case::D4(vec!(Point::new(5.0, 0.0, 0.0), Point::new(-5.0, 0.0, 0.0), Point::new(0.0, 0.0, 5.0), Point::new(0.0, 0.0, -5.0), 
        Point::new(2.0, 1.0, 3.0), Point::new(3.0, -1.0, 2.0), Point::new(-2.0, -1.0, 3.0), Point::new(-3.0, 1.0, 2.0), Point::new(2.0, -1.0, -3.0),
        Point::new(3.0, 1.0, -2.0), Point::new(-2.0, 1.0, -3.0), Point::new(-3.0, -1.0, -2.0)), 3, "D 4")]
    #[case::D5(vec!(Point::new(0.0, 1.0, 1.0), Point::new((2.0 * PI / 5.0).sin(), (2.0 * PI / 5.0).cos(), 1.0), 
        Point::new((4.0 * PI / 5.0).sin(), -(PI / 5.0).cos(), 1.0),  Point::new(-(4.0 * PI / 5.0).sin(), -(PI / 5.0).cos(), 1.0), 
        Point::new(-(2.0 * PI / 5.0).sin(), (2.0 * PI / 5.0).cos(), 1.0), Point::new(0.0, 0.0, 0.0), 
        Point::new(0.7071067811865476, -0.7071067811865475, -1.0), 
        Point::new(0.8910065241883678, 0.45399049973954675, -1.0), Point::new(-0.1564344650402309, 0.9876883405951378, -1.0),  
        Point::new(-0.9876883405951378, 0.15643446504023079, -1.0), Point::new(-0.45399049973954675, -0.8910065241883678, -1.0)), 4, "D 5")]
    #[case::D6(vec!(Point::new(1.0, 0.0, 1.0), Point::new(-1.0, 0.0, 1.0), Point::new(0.5, 3_f64.sqrt()/2.0, 1.0), Point::new(0.0, 0.0, 0.0),
        Point::new(-0.5, 3_f64.sqrt()/2.0, 1.0), Point::new(-0.5, -3_f64.sqrt()/2.0, 1.0), Point::new(0.5, -3_f64.sqrt()/2.0, 1.0), 
        Point::new(0.25881904510252074, 0.9659258262890682, -1.0), Point::new(-0.7071067811865476, 0.7071067811865475, -1.0), 
        Point::new(-0.9659258262890682, -0.2588190451025209, -1.0), Point::new(-0.258819045102521, -0.9659258262890681, -1.0), 
        Point::new(0.7071067811865475, -0.7071067811865476, -1.0), Point::new(0.9659258262890683, 0.2588190451025203, -1.0)), 2, "D 6")]
    #[case::D7(vec!(Point::new(0.0, 1.0, 1.0), Point::new((2.0 * PI / 7.0).sin(), (2.0 * PI / 7.0).cos(), 1.0), Point::new((4.0 * PI / 7.0).sin(), (4.0 * PI / 7.0).cos(), 1.0), 
        Point::new((6.0 * PI / 7.0).sin(), (6.0 * PI / 7.0).cos(), 1.0), Point::new((8.0 * PI / 7.0).sin(), (8.0 * PI / 7.0).cos(), 1.0), 
        Point::new((10.0 * PI / 7.0).sin(), (10.0 * PI / 7.0).cos(), 1.0), Point::new((12.0 * PI / 7.0).sin(), (12.0 * PI / 7.0).cos(), 1.0),
        Point::new(0.30901699437494745, -0.9510565162951536, -1.0), Point::new(0.9362348706397372, -0.35137482408134285, -1.0), 
        Point::new(0.8584487936018662, 0.5128992774059061, -1.0), Point::new(0.13423326581765566, 0.9909497617679347, -1.0), 
        Point::new(-0.6910626489868646, 0.7227948638273916, -1.0), Point::new(-0.9959742939952392, -0.08963930890343336, -1.0), 
        Point::new(-0.5508969814521026, -0.8345732537213026, -1.0), Point::new(0.0, 0.0, 0.0)), 4, "D 7")]
    #[case::D8(vec!(Point::new(1.0, 0.0, 1.0), Point::new(2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 1.0), Point::new(0.0, 1.0, 1.0), 
        Point::new(-2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 1.0), Point::new(-1.0, 0.0, 1.0), Point::new(-2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 1.0), 
        Point::new(0.0, -1.0, 1.0), Point::new(2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 1.0), Point::new(0.0, 0.0, 0.0),
        Point::new(0.07845909572784499, 0.9969173337331281, -1.0), Point::new(0.760405965600031, 0.6494480483301838, -1.0), 
        Point::new(0.9969173337331281, -0.07845909572784493, -1.0), Point::new(0.6494480483301839, -0.760405965600031, -1.0), 
        Point::new(-0.07845909572784482, -0.9969173337331281, -1.0), Point::new(-0.7604059656000306, -0.6494480483301841, -1.0), 
        Point::new(-0.9969173337331281, 0.07845909572784521, -1.0), Point::new(-0.6494480483301842, 0.7604059656000306, -1.0)), 2, "D 8")]
    #[case::C2h(vec!(Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new(2.0, 0.0, 1.0), Point::new(-2.0, 0.0, -1.0)), 2, "C 2 h")]
    #[case::C3h(vec!(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 3_f64.sqrt()/3.0, 0.0), Point::new(-0.5, -3_f64.sqrt()/6.0, 0.0), 
        Point::new(0.5, -3_f64.sqrt()/6.0, 0.0), Point::new(1.5*-0.17841104488654497, 1.5*0.5490927356975547, 0.0), 
        Point::new(1.5*-0.3863227357043043, 1.5*-0.42905486503625107, 0.0), Point::new(1.5*0.5647337805908493, 1.5*-0.12003787066130361, 0.0)), 3, "C 3 h")]
    #[case::C4h(generate_Cnh_coordinates(4.0, false), 2, "C 4 h")]
    #[case::C5h(generate_Cnh_coordinates(5.0, true), 2, "C 5 h")]
    #[case::C6h(generate_Cnh_coordinates(6.0, true), 2, "C 6 h")]
    #[case::C7h(generate_Cnh_coordinates(7.0, true), 2, "C 7 h")]
    #[case::C8h(generate_Cnh_coordinates(8.0, true), 2, "C 8 h")]
    #[case::C2v_water(vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)], 2, "C 2 v")]
    #[case::C2v_cyclohexane(vec!(Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0), Point::new(0.0, -1.0, 0.0), 
        Point::new(1.0, 1.0, 0.5), Point::new(-1.0, -1.0, 0.5), Point::new(1.0, 0.0, -1.0), Point::new(-1.0, 0.0, -1.0), 
        Point::new(0.0, 1.0, -1.0), Point::new(0.0, -1.0, -1.0), Point::new(1.5, 0.0, 0.5), Point::new(-1.5, 0.0, 0.5), 
        Point::new(0.0, 1.5, 0.5), Point::new(0.0, -1.5, 0.5), Point::new(1.0, 1.0, 1.5), Point::new(-1.0, -1.0, 1.5),
        Point::new(2.0, 2.0, 0.5), Point::new(-2.0, -2.0, 0.5)), 2, "C 2 v")]
    #[case::C3v_ammonia(vec!(Point::new(0.0, 0.0, 1.0), Point::new(0.0, 3_f64.sqrt()/3.0, 0.0), Point::new(-0.5, -3_f64.sqrt()/6.0, 0.0), 
        Point::new(0.5, -3_f64.sqrt()/6.0, 0.0)), 3, "C 3 v")]
    #[case::C3v_triflate(vec!(Point::new(0.0, 0.0, 1.0), Point::new(0.0, 3_f64.sqrt()/3.0, 1.0), Point::new(-0.5, -3_f64.sqrt()/6.0, 1.0), 
        Point::new(0.5, -3_f64.sqrt()/6.0, 1.0), Point::new(0.0, 0.0, -1.0), Point::new(0.0, -3_f64.sqrt()/3.0, -1.5), Point::new(-0.5, 3_f64.sqrt()/6.0, -1.5), 
        Point::new(0.5, 3_f64.sqrt()/6.0, -1.5)), 3, "C 3 v")]
    #[case::C3v_complex(vec!(Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new(0.0, 2.0, 0.0), 
        Point::new(0.0, -1.0, 0.0), Point::new(0.0, 0.0, 2.0), Point::new(0.0, 0.0, -1.0)), 3, "C 3 v")]
    #[case::C4v(vec!(Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0), 
        Point::new(0.0, -1.0, 0.0), Point::new(0.0, 0.0, 1.0), Point::new(0.0, 0.0, -1.0)), 3, "C 4 v")]
    #[case::C5v(generate_Cnv_coordinates(5.0), 2, "C 5 v")]
    #[case::C6v(generate_Cnv_coordinates(6.0), 2, "C 6 v")]
    #[case::C7v(generate_Cnv_coordinates(7.0), 2, "C 7 v")]
    #[case::C8v(generate_Cnv_coordinates(8.0), 2, "C 8 v")]
    #[case::S4(vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0), Point::new(-1.0, 1.0, -1.0), Point::new(-1.0, -1.0, 1.0), 
        Point::new(1.0, -1.0, -1.0), Point::new(0.5, 1.2, 0.5), Point::new(-0.5, 1.2, -0.5), Point::new(0.5, -1.2, -0.5),
        Point::new(-0.5, -1.2, 0.5), Point::new(2.0, 1.0, 1.0), Point::new(-2.0, 1.0, -1.0), Point::new(-1.0, -1.0, 2.0), 
        Point::new(1.0, -1.0, -2.0)), 4, "S 4")]
    #[case::Cs(vec!(Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new(1.5, -1.0, 0.0), Point::new(-1.5, -1.0, 0.0),
        Point::new(1.5, 1.5, 0.0), Point::new(-1.5, 2.0, 0.0)), 2, "C s")]
    #[case::Ci(vec!(Point::new(0.0, 0.0, 1.0), Point::new(0.0, 2.0*3_f64.sqrt()/3.0, 1.0), Point::new(-0.5*3.0, -3_f64.sqrt()/6.0*3.0, 1.0), 
        Point::new(0.5, -3_f64.sqrt()/6.0, 1.0), Point::new(0.0, 0.0, -1.0), Point::new(0.0, -2.0*3_f64.sqrt()/3.0, -1.0), Point::new(-0.5, 3_f64.sqrt()/6.0, -1.0), 
        Point::new(0.5*3.0, 3_f64.sqrt()/6.0*3.0, -1.0)), 3, "C i")]
    #[case::C1(vec!(Point::new(0.0, 0.0, 1.0), Point::new(0.0, 2.0*3_f64.sqrt()/3.0, 1.0), Point::new(-0.5*3.0, -3_f64.sqrt()/6.0*3.0, 1.0), 
        Point::new(0.5, -3_f64.sqrt()/6.0, 1.0), Point::new(0.0, 0.0, -1.0), Point::new(0.0, -2.0*3_f64.sqrt()/3.0, -1.0), Point::new(-0.5, 3_f64.sqrt()/6.0, -1.0), 
        Point::new(0.5, 3_f64.sqrt()/6.0, -1.0)), 3, "C 1")]
    fn test_find_point_group(#[case] coordinates: Vec<Point>, #[case] precision: u16, #[case] expected: String) {
        let result = find_point_group(coordinates, precision, 1e-10);
        assert_eq!(expected, result)
    }

    #[test]
    fn test_is_linear_true() {
        assert!(is_linear(&vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)), 1e-10))
    }

    #[test]
    fn test_is_linear_false() {
        let result = is_linear(&vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)), 1e-10);
        assert!(!result, "Expected value: false, Actual value: {}", result)
    }

    #[test]
    fn test_get_rotation_axes_precision_1() {
        let coordinates = vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0), Point::new(2.0, 0.0, 0.0)];
        let result = _get_rotation_axes_thorough(&coordinates, 1, 1e-10);

        let expected_result = vec![Line::new(coordinates[0], Vector::new(1.0, 0.0, 0.0)), 
                                   Line::new(coordinates[0], Vector::new(0.0, 1.0, 0.0)),
                                   Line::new(coordinates[1], Vector::new(-1.0, 1.0, 0.0)),
                                   Line::new(coordinates[2], Vector::new(2.0, -1.0, 0.0))];
        
        for (expected, actual) in std::iter::zip(expected_result, result) {
            assert_eq!(expected, actual)
        }
    }

    fn expected_get_rotation_axes_octagon_precision_2() -> Vec<Line> {
        let p = vec!(Point::new(1.0, 0.0, 0.0), Point::new((2.0+2_f64.sqrt())/4.0, 2_f64.sqrt()/4.0, 0.0), Point::new(0.5, 0.5, 0.0), // 0 1 2
            Point::new((2.0-2_f64.sqrt())/4.0, 2_f64.sqrt()/4.0, 0.0), Point::new(0.0, 0.0, 0.0), Point::new((2.0-2_f64.sqrt())/4.0, -2_f64.sqrt()/4.0, 0.0), // 3 4 5
            Point::new(0.5, -0.5, 0.0), Point::new((2.0+2_f64.sqrt())/4.0, -2_f64.sqrt()/4.0, 0.0), Point::new(2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 0.0), // 6 7 8
            Point::new(2_f64.sqrt()/4.0, (2.0+2_f64.sqrt())/4.0, 0.0), Point::new(0.0, 2_f64.sqrt()/2.0, 0.0), Point::new((2_f64.sqrt()-2.0)/4.0, 2_f64.sqrt()/4.0, 0.0), // 9 10 11
            Point::new(0.0, 0.0, 0.0), Point::new(2_f64.sqrt()/4.0, (2_f64.sqrt()-2.0)/4.0, 0.0), Point::new(2_f64.sqrt()/2.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0), // 12 13 14 15
            Point::new(-2_f64.sqrt()/4.0, (2.0+2_f64.sqrt())/4.0, 0.0), Point::new(-0.5, 0.5, 0.0), Point::new(-2_f64.sqrt()/4.0, (2.0-2_f64.sqrt())/4.0, 0.0), // 16 17 18
            Point::new(0.0, 0.0, 0.0), Point::new(2_f64.sqrt()/4.0, (2.0-2_f64.sqrt())/4.0, 0.0), Point::new(-2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 0.0), // 19 20 21
            Point::new((-2_f64.sqrt()-2.0)/4.0, 2_f64.sqrt()/4.0, 0.0), Point::new(-2_f64.sqrt()/2.0, 0.0, 0.0), Point::new(-2_f64.sqrt()/4.0, (-2.0+2_f64.sqrt())/4.0, 0.0),  // 22 23 24
            Point::new(0.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new((-2.0-2_f64.sqrt())/4.0, -2_f64.sqrt()/4.0, 0.0), Point::new(-0.5, -0.5, 0.0), //25 26 27 28
            Point::new((2_f64.sqrt()-2.0)/4.0, -2_f64.sqrt()/4.0, 0.0), Point::new(-2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 0.0), //29 30
            Point::new(-2_f64.sqrt()/4.0, (-2.0-2_f64.sqrt())/4.0, 0.0), Point::new(0.0, -2_f64.sqrt()/2.0, 0.0), Point::new(0.0, -1.0, 0.0),  //31 32 33
            Point::new(2_f64.sqrt()/4.0, (-2.0-2_f64.sqrt())/4.0, 0.0), Point::new(2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 0.0)); //34 35

        let mut expected = Vec::with_capacity(100);

        for (i, j) in [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 9), (0, 10), (0, 11), (0, 13), (0, 16), (0, 17), (0, 18), (0, 20),
                                     (0, 22), (0, 24), (0, 27), (0, 28), (0, 29), (0, 31), (0, 32), (0, 34), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 9),
                                     (1, 15), (1, 17), (1, 18), (1, 21), (1, 23), (1, 26), (1, 28), (1, 29), (1, 30), (1, 32), (1, 33),(1, 35), (2, 3), (2, 4), (2, 5), 
                                     (2, 6), (2, 7), (2, 11), (2, 13), (2, 17), (2, 21), (2, 22), (2, 24), (2, 26), (2, 27), (2, 29), (2, 31), (2, 33), (2, 34), (2, 35), 
                                     (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 14), (3, 15), (3, 17), (3, 24), (3, 26), (3, 27), (3, 28), (3, 30), (3, 32), (3, 33), 
                                     (3, 34), (3, 35), (4, 5), (4, 6), (4, 7), (4, 10), (5, 6), (5, 7), (5, 8), (5, 9), (5, 10), (5, 14), (5, 15), (5, 17), (5, 18), 
                                     (5, 21), (5, 22), (5, 26), (5, 28), (5, 33), (5, 35), (6, 7), (6, 8), (6, 9), (6, 11), (6, 15), (6, 16), (6, 18), (6, 20),
                                     (6, 22), (6, 26), (6, 27), (6, 28), (6, 29), (6, 30), (7, 8), (7, 10), (7, 11), (7, 15), (7, 21), (7, 23), (7, 24), (7, 26),
                                     (7, 28), (7, 30), (7, 33), (7, 34), (8, 9), (8, 10), (8, 11), (8, 13), (8, 14), (8, 16), (8, 17), (8, 18), (8, 20), (8, 22), (8, 23),
                                     (8, 24), (8, 27), (8, 29), (8, 31), (8, 32), (8, 34), (9, 10), (9, 11), (9, 13), (9, 16), (9, 21), (9, 23), (9, 24), (9, 26), (9, 27),
                                     (9, 28), (9, 30), (9, 32), (9, 33), (9, 35), (10, 11), (10, 14), (10, 16), (10, 18), (10, 20), (10, 23), (10, 26), (10, 27), 
                                     (10, 29), (10, 30), (10, 31), (11, 13), (11, 14), (11, 21), (11, 23), (11, 29), (11, 30), (11, 31), (11, 32), (11, 33), (11, 35),
                                     (13, 14), (13, 15), (13, 16), (13, 17), (13, 21), (13, 23), (13, 24), (13, 26), (13, 27), (13, 30), (13, 32), (13, 35), (14, 15),
                                     (14, 16), (14, 18), (14, 21), (14, 22), (14, 24), (14, 27), (14, 28), (14, 30), (14, 31), (14, 32), (14, 33), (15, 16), (15, 17),
                                     (15, 18), (15, 20), (15, 22), (15, 23), (15, 24), (15, 27), (15, 28), (15, 29), (15, 31), (15, 34), (16, 17), (16, 18), (16, 22),
                                     (16, 26), (16, 28), (16, 29), (16, 30), (16, 32), (16, 33), (16, 35), (17, 18), (17, 24), (17, 28), (17, 30), (17, 31), (17, 33),
                                     (17, 34), (18, 20), (18, 21), (18, 26), (18, 33), (18, 34), (18, 35), (20, 21), (20, 22), (20, 23), (20, 26), (20, 28), (20, 29),
                                     (20, 30), (20, 31), (20, 33), (21, 22), (21, 23), (21, 24), (21, 27), (21, 28), (21, 29), (21, 31), (21, 32), (21, 34), (22, 23),
                                     (22, 24), (22, 27), (22, 30), (22, 32), (22, 33), (22, 35), (23, 29), (23, 32), (23, 33), (23, 34), (23, 35), (24, 26), (24, 30),
                                     (24, 32), (24, 35), (26, 27), (26, 28), (26, 29), (26, 31), (26, 32), (26, 34), (27, 28), (27, 31), (27, 33), (27, 35), (28, 35),
                                     (29, 30), (29, 33), (30, 31), (30, 32), (30, 34), (31, 34), (31, 35), (33, 34)] {
            match Line::from_points(p[i], p[j]) {
                Some(line) => expected.push(line),
                None => ()
            }
        }

        return expected
    }

    #[rstest]
    #[ignore]
    #[case::octagon2(octagon_coordinates(), 2, expected_get_rotation_axes_octagon_precision_2())]
    fn test_get_rotation_axes(#[case] coordinates: Vec<Point>, #[case] precision: u16, #[case] expected: Vec<Line>) {
        let result = _get_rotation_axes_thorough(&coordinates, precision, 1e-10);
        assert_eq!(expected, result);
    }

    #[rstest]
    #[case(vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)), 2, 
    vec!(Point::new(0.0, 0.0, 0.0), Point::new(0.5, 0.0, 0.0), Point::new(0.0, 0.5, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.5, 0.5, 0.0), Point::new(0.0, 1.0, 0.0)))]
    #[case(vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)), 3, 
    vec!(Point::new(0.0, 0.0, 0.0), Point::new(0.5, 0.0, 0.0), Point::new(1.0/3.0, 1.0/3.0, 0.0), Point::new(0.0, 0.5, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.5, 0.5, 0.0), Point::new(0.0, 1.0, 0.0)))]
    #[case::octagon2(octagon_coordinates(), 
        2, 
        vec!(Point::new(1.0, 0.0, 0.0), Point::new((2.0+2_f64.sqrt())/4.0, 2_f64.sqrt()/4.0, 0.0), Point::new(0.5, 0.5, 0.0), 
            Point::new((2.0-2_f64.sqrt())/4.0, 2_f64.sqrt()/4.0, 0.0), Point::new(0.0, 0.0, 0.0), Point::new((2.0-2_f64.sqrt())/4.0, -2_f64.sqrt()/4.0, 0.0),
            Point::new(0.5, -0.5, 0.0), Point::new((2.0+2_f64.sqrt())/4.0, -2_f64.sqrt()/4.0, 0.0), Point::new(2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 0.0),
            Point::new(2_f64.sqrt()/4.0, (2.0+2_f64.sqrt())/4.0, 0.0), Point::new(0.0, 2_f64.sqrt()/2.0, 0.0), Point::new((2_f64.sqrt()-2.0)/4.0, 2_f64.sqrt()/4.0, 0.0),
            Point::new(0.0, 0.0, 0.0), Point::new(2_f64.sqrt()/4.0, (2_f64.sqrt()-2.0)/4.0, 0.0), Point::new(2_f64.sqrt()/2.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0), 
            Point::new(-2_f64.sqrt()/4.0, (2.0+2_f64.sqrt())/4.0, 0.0), Point::new(-0.5, 0.5, 0.0), Point::new(-2_f64.sqrt()/4.0, (2.0-2_f64.sqrt())/4.0, 0.0),
            Point::new(0.0, 0.0, 0.0), Point::new(2_f64.sqrt()/4.0, (2.0-2_f64.sqrt())/4.0, 0.0), Point::new(-2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 0.0), 
            Point::new((-2_f64.sqrt()-2.0)/4.0, 2_f64.sqrt()/4.0, 0.0), Point::new(-2_f64.sqrt()/2.0, 0.0, 0.0), Point::new(-2_f64.sqrt()/4.0, (-2.0+2_f64.sqrt())/4.0, 0.0), 
            Point::new(0.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), Point::new((-2.0-2_f64.sqrt())/4.0, -2_f64.sqrt()/4.0, 0.0), Point::new(-0.5, -0.5, 0.0), 
            Point::new((2_f64.sqrt()-2.0)/4.0, -2_f64.sqrt()/4.0, 0.0), Point::new(-2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 0.0), 
            Point::new(-2_f64.sqrt()/4.0, (-2.0-2_f64.sqrt())/4.0, 0.0), Point::new(0.0, -2_f64.sqrt()/2.0, 0.0), Point::new(0.0, -1.0, 0.0), 
            Point::new(2_f64.sqrt()/4.0, (-2.0-2_f64.sqrt())/4.0, 0.0), Point::new(2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 0.0)))]

    #[case::octagon3(octagon_coordinates(), 
        3, 
        vec!(Point::new(1.0, 0.0, 0.0), Point::new((2.0+2_f64.sqrt())/4.0, 2_f64.sqrt()/4.0, 0.0), Point::new((2.0+2_f64.sqrt())/6.0, (2.0+2_f64.sqrt())/6.0, 0.0), 
            Point::new(1.0/3.0, 2_f64.sqrt()/3.0, 0.0), Point::new(2_f64.sqrt()/6.0, 2_f64.sqrt()/6.0, 0.0), Point::new(1.0/3.0, 0.0, 0.0), 
            Point::new((2.0+2_f64.sqrt())/6.0, (-2.0+2_f64.sqrt())/6.0, 0.0), Point::new((1.0+2_f64.sqrt())/3.0, 0.0, 0.0), Point::new(0.5, 0.5, 0.0), 
            Point::new((2.0-2_f64.sqrt())/6.0, (2.0+2_f64.sqrt())/6.0, 0.0), Point::new(0.0, 1.0/3.0, 0.0), Point::new((2.0-2_f64.sqrt())/6.0, (2.0-2_f64.sqrt())/6.0, 0.0),
            Point::new(1.0/3.0, 0.0, 0.0), Point::new((2.0+2_f64.sqrt())/6.0, (2.0-2_f64.sqrt())/6.0, 0.0), Point::new((2.0-2_f64.sqrt())/4.0, 2_f64.sqrt()/4.0, 0.0), 
            Point::new(-2_f64.sqrt()/6.0, 2_f64.sqrt()/6.0, 0.0), Point::new((1.0-2_f64.sqrt())/3.0, 0.0, 0.0), Point::new((2.0-2_f64.sqrt())/6.0, (-2.0+2_f64.sqrt())/6.0, 0.0),
            Point::new(1.0/3.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0), Point::new(-2_f64.sqrt()/6.0, -2_f64.sqrt()/6.0, 0.0), Point::new(0.0, -1.0/3.0, 0.0), 
            Point::new(2_f64.sqrt()/6.0, -2_f64.sqrt()/6.0, 0.0), Point::new((2.0-2_f64.sqrt())/4.0, -2_f64.sqrt()/4.0, 0.0), Point::new((2.0-2_f64.sqrt())/6.0,(-2.0-2_f64.sqrt())/6.0, 0.0 ),
            Point::new(1.0/3.0, -2_f64.sqrt()/3.0, 0.0), Point::new(0.5, -0.5, 0.0), Point::new((2.0+2_f64.sqrt())/6.0, (-2.0-2_f64.sqrt())/6.0, 0.0),
            Point::new((2.0+2_f64.sqrt())/4.0, -2_f64.sqrt()/4.0, 0.0), Point::new(2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 0.0), Point::new(2_f64.sqrt()/4.0, (2.0+2_f64.sqrt())/4.0, 0.0), 
            Point::new(0.0, (1.0+2_f64.sqrt())/3.0, 0.0), Point::new((-2.0+2_f64.sqrt())/6.0, (2.0+2_f64.sqrt())/6.0, 0.0), Point::new(0.0, 1.0/3.0, 0.0),
            Point::new(2_f64.sqrt()/6.0, 2_f64.sqrt()/6.0, 0.0), Point::new(2_f64.sqrt()/3.0, 1.0/3.0, 0.0), Point::new(0.0, 2_f64.sqrt()/2.0, 0.0), 
            Point::new(-1.0/3.0, 2_f64.sqrt()/3.0, 0.0), Point::new(-2_f64.sqrt()/6.0, 2_f64.sqrt()/6.0, 0.0), Point::new(0.0, (-1.0+2_f64.sqrt())/3.0, 0.0),
            Point::new(2_f64.sqrt()/6.0, 2_f64.sqrt()/6.0, 0.0), Point::new((2_f64.sqrt()-2.0)/4.0, 2_f64.sqrt()/4.0, 0.0), Point::new(-1.0/3.0, 0.0, 0.0), 
            Point::new((-2.0+2_f64.sqrt())/6.0, (-2.0+2_f64.sqrt())/6.0, 0.0), Point::new((-1.0+2_f64.sqrt())/3.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0), 
            Point::new(0.0, -1.0/3.0, 0.0), Point::new(2_f64.sqrt()/6.0, -2_f64.sqrt()/6.0, 0.0), Point::new(2_f64.sqrt()/4.0, (2_f64.sqrt()-2.0)/4.0, 0.0), 
            Point::new(2_f64.sqrt()/3.0, -1.0/3.0, 0.0), Point::new(2_f64.sqrt()/2.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0), Point::new(-2_f64.sqrt()/4.0, (2.0+2_f64.sqrt())/4.0, 0.0), 
            Point::new((-2.0-2_f64.sqrt())/6.0, (2.0+2_f64.sqrt())/6.0, 0.0), Point::new(-2_f64.sqrt()/3.0, 1.0/3.0, 0.0), Point::new(-2_f64.sqrt()/6.0, 2_f64.sqrt()/6.0, 0.0),
            Point::new(0.0, 1.0/3.0, 0.0), Point::new(-0.5, 0.5, 0.0), Point::new((-2.0-2_f64.sqrt())/6.0, (2.0-2_f64.sqrt())/6.0, 0.0), Point::new(-1.0/3.0, 0.0, 0.0),
            Point::new((-2.0+2_f64.sqrt())/6.0, (2.0-2_f64.sqrt())/6.0, 0.0), Point::new(-2_f64.sqrt()/4.0, (2.0-2_f64.sqrt())/4.0, 0.0),
            Point::new(-2_f64.sqrt()/6.0, -2_f64.sqrt()/6.0, 0.0), Point::new(0.0, (1.0-2_f64.sqrt())/3.0, 0.0), Point::new(0.0, 0.0, 0.0), 
            Point::new(2_f64.sqrt()/6.0, -2_f64.sqrt()/6.0, 0.0), Point::new(2_f64.sqrt()/4.0, (2.0-2_f64.sqrt())/4.0, 0.0), Point::new(-2_f64.sqrt()/2.0, 2_f64.sqrt()/2.0, 0.0), 
            Point::new((-2_f64.sqrt()-2.0)/4.0, 2_f64.sqrt()/4.0, 0.0), Point::new((-1.0-2_f64.sqrt())/3.0, 0.0, 0.0), Point::new((-2.0-2_f64.sqrt())/6.0, (-2.0+2_f64.sqrt())/6.0, 0.0),
            Point::new(-1.0/3.0, 0.0, 0.0), Point::new(-2_f64.sqrt()/2.0, 0.0, 0.0), Point::new(-2_f64.sqrt()/3.0, -1.0/3.0, 0.0), Point::new(-2_f64.sqrt()/6.0, -2_f64.sqrt()/6.0, 0.0),
            Point::new(-2_f64.sqrt()/4.0, (-2.0+2_f64.sqrt())/4.0, 0.0), Point::new(0.0, -1.0/3.0, 0.0), Point::new(0.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0), 
            Point::new((-2.0-2_f64.sqrt())/4.0, -2_f64.sqrt()/4.0, 0.0), Point::new((-2.0-2_f64.sqrt())/6.0, (-2.0-2_f64.sqrt())/6.0, 0.0), Point::new(-1.0/3.0, -2_f64.sqrt()/3.0, 0.0),
            Point::new(-0.5, -0.5, 0.0), Point::new((-2.0+2_f64.sqrt())/6.0, (-2.0-2_f64.sqrt())/6.0, 0.0), Point::new((2_f64.sqrt()-2.0)/4.0, -2_f64.sqrt()/4.0, 0.0), 
            Point::new(-2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 0.0), Point::new(-2_f64.sqrt()/4.0, (-2.0-2_f64.sqrt())/4.0, 0.0), Point::new(0.0, (-1.0-2_f64.sqrt())/3.0, 0.0),
            Point::new(0.0, -2_f64.sqrt()/2.0, 0.0), Point::new(0.0, -1.0, 0.0), Point::new(2_f64.sqrt()/4.0, (-2.0-2_f64.sqrt())/4.0, 0.0), Point::new(2_f64.sqrt()/2.0, -2_f64.sqrt()/2.0, 0.0)))]
    fn test_create_pseudo_points(#[case] coordinates: Vec<Point>, #[case] precision: u16, #[case] expected: Vec<Point>) {
        let result = create_pseudo_points(&coordinates, precision, 1e-10);
        assert_eq!(expected, result);
    }

    #[test]
    fn test_recompute_baseline_rotation_8_to_4() {
        let coordinates = vec![Point::new(0.0, 1.0, 0.0), Point::new(0.0, 0.0, 1.0), Point::new(0.0, 0.0, -1.0)];
        let axis = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
        let result = recompute_baseline_rotation(&coordinates, &Point::new(0.0, 1.0, 0.0), &axis, 8.0, 1e-10);

        assert_eq!(4.0, result)
    }

    #[test]
    fn test_recompute_baseline_rotation_n_to_2() {
        let coordinates = vec![Point::new(0.0, 1.0, 0.0), Point::new(0.0, -1.0, 0.0)];
        let axis = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));

        for i in [8.0, 4.0, 6.0] {
            let result = recompute_baseline_rotation(&coordinates, &Point::new(0.0, 1.0, 0.0), &axis, i, 1e-10);

            assert_eq!(2.0, result)
        }
    }
    
    #[test]
    fn test_recompute_baseline_rotation_n_to_1() {
        let coordinates = vec![Point::new(0.0, 1.0, 0.0), Point::new(0.0, 0.0, 3.0), Point::new(0.0, 0.0, -3.0)];
        let axis = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));

        for i in 2..=8 {
            let result = recompute_baseline_rotation(&coordinates, &Point::new(0.0, 1.0, 0.0), &axis, i as f64, 1e-10);

            assert_eq!(1.0, result)
        }
    }

    #[test]
    fn test_recompute_baseline_rotation_6_to_3() {
        let coordinates = vec![Point::new(0.0, 1.0, 0.0), Point::new(0.0, -1.0, 0.0), Point::new(1.0, -0.5, 0.8660254037844388)];
        let axis = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
        let result = recompute_baseline_rotation(&coordinates, &Point::new(1.0, 1.0, 0.0), &axis, 6.0, 1e-10);

        assert_eq!(3.0, result)
    }

    #[test]
    fn test_is_rotation_viable_true() {
        let coordinates = vec![Point::new(1.0, 1.0, 1.0), Point::new(1.0, -1.0, -1.0)];
        let point = Point::new(1.0, 1.0, 1.0);
        let axis = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
        let result = is_rotation_viable(&coordinates, &point, &axis, 2.0, 1e-10);

        assert!(result, "rotated point = {:?}", axis.rotate_point(point, 2.0))
    }

    #[test]
    fn test_is_rotation_viable_false() {
        let coordinates = vec![Point::new(1.0, 1.0, 1.0), Point::new(1.0, -1.0, 1.0)];
        let point = Point::new(1.0, 1.0, 1.0);
        let axis = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
        let result = is_rotation_viable(&coordinates, &point, &axis, 2.0, 1e-10);

        assert!(!result, "rotated point = {:?}", axis.rotate_point(point, 2.0))
    }

    #[test]
    fn test_is_equivalent_point() {
        let coordinates = vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(10.0, 5.0, -5.0), Point::new(0.0, 0.0, 5.5)];

        assert!(is_equivalent_point(&coordinates, Point::new(0.0, 0.0, 0.0), 1e-10));
        assert!(is_equivalent_point(&coordinates, Point::new(0.0, 0.0, 5.5), 1e-10));
        assert!(!is_equivalent_point(&coordinates, Point::new(0.0, 0.0, 6.5), 1e-10))
    }

    #[test]
    fn test_combinations() {
        assert_eq!(10, combinations(5, 2));
        assert_eq!(15, combinations(6, 2))
    }
}
