use crate::geometry::*;


/// Computes the distance between a Point and a Plane for 
/// a set of Points. The returned distance is positive
/// when a Point is on the same side of the Plane as the
/// normal Vector (i.e. when the Vector between the Plane
/// and the Point is a positive k-multiple of the normal
/// Vector) and vice versa.
/// 
/// # Example
/// 
/// ```
/// use point_group::geometry::*;
/// use point_group::utils::*;
/// 
/// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
/// 
/// let coordinates = vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0),
///     Point::new(0.0, 0.0, 5.0), Point::new(5.0, 0.0, 0.1),
///     Point::new(0.0, 5.0, -0.0001), Point::new(-3.0, -4.0, -5.0));
/// 
/// assert_eq!(vec![0.0, 0.0, 5.0, 0.1, -0.0001, -5.0], 
///     deviation_from_plane(&coordinates, plane))
/// ```
pub fn deviation_from_plane(coordinates: &Vec<Point>, plane: Plane) -> Vec<f64> {
    let mut deviation = Vec::with_capacity(coordinates.len());
    let magnitude = plane.normal.magnitude();

    for point in coordinates {
        deviation.push(Vector::from_two_points(plane.point, *point).dot_product(plane.normal) / magnitude)
    }

    return deviation
}


/// Computes the deviation of a set of Points from a
/// rotation symmetry. In other words, it calculates
/// how far a given rotation is from being perfectly
/// symmetrical for a set of Points.
/// 
/// This is done by rotating all the Points around 
/// the given Line, and then calculating how far
/// each Point is from its closest equivalent in
/// the original set.
/// 
/// # Example
/// 
/// ```
/// use point_group::geometry::*;
/// use point_group::utils::*;
/// 
/// let axis = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
/// 
/// let coordinates = vec!(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 1.0),
///     Point::new(0.001, 0.0, 5.0), Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0),
///     Point::new(0.5, 0.0, -1.0), Point::new(-0.5, 0.0, -1.0), Point::new(0.0, 0.5, -1.0),
///     Point::new(0.0, -0.5, -1.0), Point::new(5.0, 0.0, 1.0), Point::new(-5.5, 0.0, 1.0));
/// 
/// let result = Vec::from(deviation_from_rotation_symmetry(&coordinates, axis, 2.0));
/// 
/// fn compare(arr1: Vec<f64>, arr2: Vec<f64>) -> Vec<bool> {
///     arr1.iter().zip(arr2).map(|(x, y)| (x - y).abs() < 1e-10).collect()
/// }
/// 
/// assert_eq!(vec![true; 11], compare(vec![0.0, 0.0, 0.002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5], result));
/// ```
pub fn deviation_from_rotation_symmetry(coordinates: &Vec<Point>, axis: Line, n: f64) -> Vec<f64> {
    let mut deviation = Vec::with_capacity(coordinates.len());

    for point in coordinates {
        let rotated = point.rotate_around_axis(axis, n);

        deviation.push(distance_from_closest_point(&coordinates, rotated))
    }

    return deviation
}


/// Finds the point from the provided set Point that is
/// the closest to the provided Point and computes their distance.
/// 
/// # Example
/// 
/// ```
/// use point_group::geometry::*;
/// use point_group::utils::*;
/// 
/// let coordinates = vec!(Point::new(0.0, 0.0, 0.0), Point::new(10.0, 0.0, 0.0),
///     Point::new(0.0, -5.0, 0.0), Point::new(1.0, 0.0, 1.0));
/// 
/// assert!((0.0 - distance_from_closest_point(&coordinates, Point::new(0.0, 0.0, 0.0)).abs() < 1e-10));
/// assert!((5.0 - distance_from_closest_point(&coordinates, Point::new(15.0, 0.0, 0.0)).abs() < 1e-10));
/// assert!((2.5 - distance_from_closest_point(&coordinates, Point::new(0.0, -2.5, 0.0)).abs() < 1e-10));
/// ```
pub fn distance_from_closest_point(coordinates: &Vec<Point>, point: Point) -> f64 {
    let mut diff = 100.0;
    for original in coordinates {
        let distance = Vector::from_two_points(point, *original).magnitude();

        if distance < diff {
            diff = distance
        }
    }

    return diff
}



/// Finds the Point from the provided set that is the closest
/// to the provided Point and computes their distance. Both
/// the distance and the index of the closest point are returned.
/// 
/// # Example
/// 
/// ```
/// use point_group::geometry::*;
/// use point_group::utils::*;
/// 
/// let coordinates = vec!(Point::new(0.0, 0.0, 0.0), Point::new(10.0, 0.0, 0.0),
///     Point::new(0.0, -5.0, 0.0), Point::new(1.0, 0.0, 1.0));
/// 
/// // This point is identical to the first point
/// let result1 = distance_from_closest_point_pairs(&coordinates, Point::new(0.0, 0.0, 0.0));
/// assert!((0.0 - result1.0).abs() < 1e-10);
/// assert_eq!(0, result1.1);
/// 
/// // This point is closest to the point at index 2, displaced by 5.0
/// let result2 = distance_from_closest_point_pairs(&coordinates, Point::new(15.0, 0.0, 0.0));
/// assert!((5.0 - result2.0).abs() < 1e-10);
/// assert_eq!(1, result2.1);
/// 
/// // This point is equally distanat from points at indices 0 and 2, displaced by 2.5
/// let result3 = distance_from_closest_point_pairs(&coordinates, Point::new(0.0, -2.5, 0.0));
/// assert!((2.5 - result3.0).abs() < 1e-10);
/// assert_eq!(0, result3.1)
/// ```
pub fn distance_from_closest_point_pairs(coordinates: &Vec<Point>, point: Point) -> (f64, usize) {
    let mut diff = 100.0;
    let mut closest = 0;

    for (i, original) in coordinates.iter().enumerate() {
        let distance = Vector::from_two_points(point, *original).magnitude();

        if distance < diff {
            diff = distance;
            closest = i
        }
    }

    return (diff, closest)
}


pub fn deviation_from_reflection_symmetry(coordinates: &Vec<Point>, plane: Plane) -> Vec<f64> {
    let mut deviation = Vec::with_capacity(coordinates.len());

    for point in coordinates {
        println!("{} {:?} -> {:?}", point.is_approx_on_plane(plane, 0.1), point, point.reflect(plane));
        deviation.push(distance_from_closest_point(&coordinates, point.reflect(plane)))
    }

    return deviation
}



pub fn deviation_from_reflection_symmetry_pairs(coordinates: &Vec<Point>, plane: Plane) -> (Vec<f64>, Vec<(usize, usize)>) {
    let mut deviation = Vec::with_capacity(coordinates.len());
    let mut pairs = Vec::with_capacity(coordinates.len());

    for (i, point) in coordinates.iter().enumerate() {
        let (diff, closest) = distance_from_closest_point_pairs(&coordinates, point.reflect(plane));
        deviation.push(diff);
        pairs.push((i, closest));
        println!("{:?} {:?} -> {:?} ~ {:?}", (i, closest), point, point.reflect(plane), coordinates[closest]);
    }

    return (deviation, pairs)
}



pub fn is_rotation_symmetrical(coordinates: &Vec<Point>, axis: Line, n: f64, tolerance: f64) -> bool {
    'outer: for point in coordinates {
        let rotated = axis.rotate_point(*point, n);

        for coord in coordinates {
            if coord.approx_eq(&rotated, tolerance) {
                continue 'outer
            }
        }

        return false
    }

    return true
}



/// Converts a vec of Points into a 2D vec of f64s,
/// containing the coordinate information. The shape of
/// the vec is 3 × n, where n is the number of points.
/// 
/// # Example
/// 
/// ```
/// use point_group::geometry::*;
/// use point_group::utils::*;
/// 
/// let coordinates = vec!(Point::new(0.0, 0.0, 0.0), Point::new(10.0, 0.0, 0.0),
///     Point::new(0.0, -5.0, 0.0), Point::new(1.0, 0.0, 1.0));
/// 
/// let expected = vec!(vec!(0.0, 0.0, 0.0), vec!(10.0, 0.0, 0.0),
///     vec!(0.0, -5.0, 0.0), vec!(1.0, 0.0, 1.0));
/// 
/// let result = points_into_coordinates(coordinates);
/// 
/// assert_eq!(3, result.len());
/// assert_eq!(4, result[0].len())
/// ```
pub fn points_into_coordinates(points: Vec<Point>) -> Vec<Vec<f64>> {
    let len = points.len();

    let mut output = vec![Vec::with_capacity(len); 3];

    for point in points {
        output[0].push(point.x);
        output[1].push(point.y);
        output[2].push(point.z)
    }

    return output
}



/// Computes the distance between the two closest Points in
/// a set of Points.
/// 
/// # Example
/// 
/// ```
/// use point_group::geometry::*;
/// use point_group::utils::*;
/// 
/// let coordinates = vec!(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), 
///     Point::new(0.0, 1.0, 0.0), Point::new(5.0, -5.0, 0.0), Point::new(-5.0, 5.0, 0.0));
/// 
/// assert!((1.0 - shortest_atom_distance(&coordinates)).abs() < 1e-10);
/// ```
pub fn shortest_atom_distance(coordinates: &Vec<Point>) -> f64 {
    let mut shortest = 100.0;
    for (i, point1) in coordinates.iter().enumerate() {
        for point2 in coordinates[i+1..].iter() {
            let distance = Vector::from_two_points(*point1, *point2).magnitude();
            if distance < shortest  {
                shortest = distance
            }
        }
    }

    return shortest
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deviation_from_plane() {
        let plane = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
        let coords = vec!(Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0), Point::new(2.0, 0.0, 5.0));

        assert_eq!(vec!(0.0, 1.0, 0.0), deviation_from_plane(&coords, plane));

        let p = Plane::new(Point::new(4.623141109, 9.024727586, 0.474048652), 
            Vector::new(-0.1, 0.99, 1.0)).unwrap();
        let c = vec!(Point::new(4.623141109, 9.024727586, 0.474048652));

        assert_eq!(vec!(0.0), deviation_from_plane(&c, p))
    }

    #[test]
    fn test_deviation_from_reflection_symmetry() {
        let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();

        let coords = vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(0.0, 0.0, 1.0), 
            Point::new(0.0, 0.0, -1.0), Point::new(0.0, 0.0, 2.0), Point::new(0.0, 0.0, 10.0), Point::new(0.0, 0.0, -4.0),
            Point::new(0.0, 0.0, 0.1), Point::new(0.0, 0.0, 0.0001), Point::new(5.0, 5.0, 5.0), Point::new(5.0, 5.0, -4.0)];

        let expected = vec![0.0, 0.0, 0.0, 0.0, 1.0, 6.0, 2.0, 0.1, 0.0001, 1.0, 1.0];

        assert_eq!(expected, deviation_from_reflection_symmetry(&coords, plane))
    }
}