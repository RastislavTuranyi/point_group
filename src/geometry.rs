//! 



pub trait ApproxEq {
    fn approx_eq(&self, other: &Self, tolerance: f64) -> bool;
}


pub trait VectorAlgebra<T> {
    type Output;

    fn cross_product(self, other: T) -> Self::Output;
}

fn _cross_product(rx: f64, ry: f64, rz: f64, lx: f64, ly: f64, lz: f64) -> (f64, f64, f64) {
    let x = ry * lz - rz * ly;
    let y = rz * lx - rx * lz;
    let z = rx * ly - ry * lx;

    return (x, y, z)
}


/// A point in 3D cartesian space, defined by its x, y, z cartesian coordinates. Equivalent to an (x, y, z) vector
/// that has its starting point in origin.
#[derive(Debug, Copy, Clone)]
pub struct Point {
    /// The x-coordinate of a Point.
    pub x: f64,
    /// The y-coordinate of a Point.
    pub y: f64,
    /// The z-coordinate of a Point.
    pub z: f64
}

impl std::ops::Add<Point> for Point {
    type Output = Point;

    fn add(self, rhs: Point) -> Self::Output {
        Point{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl<'a> std::ops::Add<&'a Point> for Point {
    type Output = Point;

    fn add(self, rhs: &Point) -> Self::Output {
        Point{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl std::ops::Add<Vector> for Point {
    type Output = Point;

    fn add(self, rhs: Vector) -> Self::Output {
        Point{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl<'a> std::ops::Add<&'a Vector> for Point {
    type Output = Point;

    fn add(self, rhs: &Vector) -> Self::Output {
        Point{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl std::ops::AddAssign<Point> for Point {
    fn add_assign(&mut self, rhs: Point) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<'a> std::ops::AddAssign<&'a Point> for Point {
    fn add_assign(&mut self, rhs: &Point) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl ApproxEq for Point {
    fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        return (self.x - other.x).abs() <= tolerance && (self.y - other.y).abs() <= tolerance && 
            (self.z - other.z).abs() <= tolerance
    }
}

impl std::ops::Div<f64> for Point {
    type Output = Point;

    fn div(self, rhs: f64) -> Self::Output {
        Point{x: self.x / rhs, y: self.y / rhs, z: self.z / rhs}
    }
}

impl std::ops::Mul<f64> for Point {
    type Output = Point;

    fn mul(self, rhs: f64) -> Self::Output {
        Point{x: self.x * rhs, y: self.y * rhs, z: self.z * rhs}
    }
}

impl std::ops::Sub<Point> for Point {
    type Output = Point;

    fn sub(self, rhs: Point) -> Self::Output {
        Point{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl<'a> std::ops::Sub<&'a Point> for Point {
    type Output = Point;

    fn sub(self, rhs: &Point) -> Self::Output {
        Point{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl<'a> std::ops::Sub<Point> for &'a Point {
    type Output = Point;

    fn sub(self, rhs: Point) -> Self::Output {
        Point{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl std::iter::Sum<Point> for Point {
    fn sum<I: Iterator<Item = Point>>(iter: I) -> Self {
        iter.fold(Point{x: 0.0, y: 0.0, z: 0.0}, |acc, x| acc + x)
    }
}

impl<'a> std::iter::Sum<&'a Point> for Point {
    fn sum<I: Iterator<Item = &'a Point>>(iter: I) -> Self {
        iter.fold(Point{x: 0.0, y: 0.0, z: 0.0}, |acc, x| acc + *x)
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        return (self.x - other.x).abs() <= 1e-10 && (self.y - other.y).abs() <= 1e-10 && (self.z - other.z).abs() <= 1e-10
    }
}

impl VectorAlgebra<Point> for Point {
    type Output = Point;

    fn cross_product(self, other: Point) -> Self::Output {
        let (x, y, z) = _cross_product(self.x, self.y, self.z, other.x, other.y, other.z);

        return Point{x, y, z}
    }
}

impl VectorAlgebra<Vector> for Point {
    type Output = Vector;

    fn cross_product(self, other: Vector) -> Self::Output {
        let (x, y, z) = _cross_product(self.x, self.y, self.z, other.x, other.y, other.z);

        return Vector{x, y, z}
    }
}

impl Point {
    /// Creates a new Point from its cartesian coordinates.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// assert_eq!(Point{x: 0.0, y: 5.0, z: -9.8}, Point::new(0.0, 5.0, -9.8))
    /// ```
    pub fn new(x: f64, y: f64, z: f64) -> Point {
        Point{x, y, z}
    }

    /// Creates a new Point from a 3-array that contains its cartesian coordinates in order of x, y, and z.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// assert_eq!(Point{x: 0.0, y: 5.0, z: -9.8}, Point::from_array([0.0, 5.0, -9.8]))
    /// ```
    pub fn from_array(arr: [f64; 3]) -> Point {
        Point{x: arr[0], y: arr[1], z: arr[2]}
    }

    /// Creates a new Point from a Vector, assuming that the start point of the vector is [0, 0, 0].
    /// 
    /// # Example 
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// assert_eq!(Point{x: 0.0, y: 5.0, z: -9.8}, Point::from_vector(Vector{x: 0.0, y: 5.0, z: -9.8}))
    /// ```
    pub fn from_vector(vector: Vector) -> Point {
        Point{x: vector.x, y: vector.y, z: vector.z}
    }

    /// Creates a new Vector from this Point, assuming that the start point of the vector is [0, 0, 0].
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point = Point::new(0.0, 5.0, -9.8);
    /// 
    /// assert_eq!(Vector{x: 0.0, y: 5.0, z: -9.8}, point.into_vector())
    /// ```
    pub fn into_vector(&self) -> Vector {
        Vector { x: self.x, y: self.y, z: self.z }
    }

    /// Returns the dot product between this Vector and another Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector1 = Vector::new(1.0, 0.0, 0.0);
    /// let vector2 = Vector::new(0.0, 1.0, 0.0);
    /// let vector3 = Vector::new(0.5, -4.1, 2.0);
    /// 
    /// assert!(vector1.dot_product(vector2).abs() < 1e-10);
    /// assert!(vector2.dot_product(vector1).abs() < 1e-10);
    /// 
    /// assert!((0.5 - vector1.dot_product(vector3).abs()) < 1e-10);
    /// assert!((-4.1 - (-1.0 * vector2.dot_product(vector3)).abs()) < 1e-10);
    /// ```
    pub fn dot_product(&self, other: Vector) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Computes the distance of this point from a Line. The following
    /// equation is used:
    ///     |v × (x-y)| / |v|
    /// where v is the directional vector of the Line, x is this Point,
    /// and y is a point on the line.
    /// 
    /// # Example
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert_eq!(0.0, Point::new(0.5, 0.0, 0.0).distance_from_line(line));
    /// assert_eq!(1.0, Point::new(0.0, 1.0, 0.0).distance_from_line(line));
    /// ```
    pub fn distance_from_line(self, line: Line) -> f64 {
        line.vector.cross_product(Vector::from_two_points(self, line.point)).magnitude() / line.vector.magnitude()
    }

    pub fn distance_from_plane(self, plane: Plane) -> f64 {
        Vector::from_two_points(plane.point, self).dot_product(plane.normal) / plane.normal.magnitude()
    }

    pub fn displacement_from_plane(self, plane: Plane) -> Vector {
        let distance = self.distance_from_plane(plane);
        return plane.normal * distance * -1.0 / plane.normal.magnitude()
    }

    /// Determines whether this Point lies on a given Line. This is a wrapper around
    /// the is_approx_on_line() method but with the tolerance set to 1e-10.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert!(Point::new(0.0, 0.0, 0.0).is_on_line(line));
    /// assert!(Point::new(5.9, 0.0, 0.0).is_on_line(line));
    /// assert!(Point::new(-0.2, 0.0, 0.0).is_on_line(line));
    /// assert!(!Point::new(0.0, 0.5, 0.0).is_on_line(line))
    /// ```
    pub fn is_on_line(&self, line: Line) -> bool {
        return self.is_approx_on_line(line, 1e-10)
    }

    /// Determines whether this Point lies on a given Line, taking into account the
    /// provided tolerance.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert!(Point::new(0.0, 0.0, 0.0).is_approx_on_line(line, 1e-10));
    /// assert!(Point::new(5.9, 0.0, 0.0).is_approx_on_line(line, 1e-5));
    /// assert!(Point::new(-0.2, 0.0, 0.0).is_approx_on_line(line, 1e-1));
    /// assert!(!Point::new(0.0, 0.5, 0.0).is_approx_on_line(line, 1e-10));
    /// 
    /// assert!(Point::new(1.0, 0.00001, 0.0).is_approx_on_line(line, 1e-4));
    /// assert!(!Point::new(1.0, 0.00001, 0.0).is_approx_on_line(line, 1e-6));
    /// ```
    pub fn is_approx_on_line(&self, line: Line, tolerance: f64) -> bool {
        let vector = Vector::from_two_points(*self, line.point);

        // If the point is the same as the first point defining the line, they are identical, otherwise check if the vectors are k multiples
        if vector.x.abs() < 1e-10 && vector.y.abs() < 1e-10 && vector.z.abs() < 1e-10 {
            return true
        } else {
            return vector.is_approx_k_multiple(line.vector, tolerance)
        }
    }

    /// Determines whether this Point lies on a given Plane. This is a
    /// wrapper around the is_approx_on_plane() method but with the
    /// tolerance set to 1e-10.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
    /// 
    /// assert!(Point::new(0.0, 0.0, 0.0).is_on_plane(plane));
    /// assert!(Point::new(9.5, 0.0, 0.0).is_on_plane(plane));
    /// assert!(Point::new(0.0, 0.0, -1.8).is_on_plane(plane));
    /// assert!(Point::new(-11.0, 0.0, 5.6).is_on_plane(plane));
    /// 
    /// assert!(!Point::new(0.0, 1.4, 0.0).is_on_plane(plane))
    /// ```
    pub fn is_on_plane(self, plane: Plane) -> bool {
        return self.is_approx_on_plane(plane, 1e-10)
    }

    /// Determines whether this Point lies on a given Plane, taking into
    /// account the provided tolerance.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
    /// 
    /// assert!(Point::new(0.0, 0.0, 0.0).is_approx_on_plane(plane, 1e-10));
    /// assert!(Point::new(9.5, 0.0, 0.0).is_approx_on_plane(plane, 1e-5));
    /// assert!(Point::new(0.0, 0.0, -1.8).is_approx_on_plane(plane, 1e-1));
    /// assert!(Point::new(-11.0, 0.0, 5.6).is_approx_on_plane(plane, 1e-13));
    /// 
    /// assert!(!Point::new(0.0, 1.4, 0.0).is_approx_on_plane(plane, 1e-10));
    /// assert!(Point::new(0.0, 1.4, 0.0).is_approx_on_plane(plane, 1e2));
    /// ```
    pub fn is_approx_on_plane(self, plane: Plane, tolerance: f64) -> bool {
        self.distance_from_plane(plane).abs() < tolerance
    }

    /// Performs an improper rotation, i.e. rotation followd by reflection, of this Point 
    /// using a rotation Quaternion and a provided mirror Plane.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let quaternion = Quaternion::new(0.0, Vector::new(1.0, 0.0, 0.0));
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(Point::new(1.0, 1.0, -1.0), Point::new(1.0, 1.0, 1.0).improper_rotate(quaternion, plane))
    /// ```
    pub fn improper_rotate(self, quaternion: Quaternion, plane: Plane) -> Point {
        self.rotate(quaternion).reflect(plane)
    }

    /// Performs an improper rotation Sn, i.e. rotation followd by reflection, of this Point 
    /// around a given Vector and a given Plane. The rotation is over 360/n degrees.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector = Vector::new(1.0, 0.0, 0.0);
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(Point::new(1.0, 1.0, -1.0), Point::new(1.0, 1.0, 1.0).improper_rotate_around_vector(vector, 2.0, plane))
    /// ```
    pub fn improper_rotate_around_vector(self, vector: Vector, n: f64, plane: Plane) -> Point {
        self.rotate_around_vector(vector, n).reflect(plane)
    }

    /// Performs an improper rotation Sn, i.e. rotation followd by reflection, of this Point 
    /// around a given Line and a given Plane. The rotation is over 360/n degrees. Please note
    /// that since the rotation occurs around an axis, the position of the Line matters.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let axis = Line::new(Point::new(0.0, 0.0, 1.0), Vector::new(1.0, 0.0, 0.0));
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Point::new(1.0, 1.0, 1.0), point.improper_rotate_around_axis(axis, 2.0, plane))
    /// ```
    pub fn improper_rotate_around_axis(self, axis: Line, n: f64, plane: Plane) -> Point {
        self.rotate_around_axis(axis, n).reflect(plane)
    }

    /// Performs an inversion of a given Point around this Point (i.e. this Point is the centre of inversion).
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let centre_of_inversion = Point::new(0.0, 0.0, 0.0);
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Point::new(-1.0, -1.0, -1.0), centre_of_inversion.invert(point))
    /// ```
    pub fn invert(self, point: Point) -> Point {
        (Vector::from_two_points(point, self) * 2.0).end_point(point)
    }

    /// Performs an inversion of this Point around a given Point (i.e. the provided Point is the centre of inversion).
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let centre_of_inversion = Point::new(0.0, 0.0, 0.0);
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Point::new(-1.0, -1.0, -1.0), point.invert_around(centre_of_inversion))
    /// ```
    pub fn invert_around(self, centre_of_inversion: Point) -> Point {
        (Vector::from_two_points(self, centre_of_inversion) * 2.0).end_point(self)
    }

    /// Performs a reflection of this Point through a given mirror Plane.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(Point::new(1.0, -1.0, 1.0), point.reflect(plane));
    /// assert_eq!(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0).reflect(plane))
    /// ```
    pub fn reflect(&self, plane: Plane) -> Point {
        // First compute the plane equation ax + by + cz + d = 0 where a, b, and c are the components of the normal
        // and d = ax0 + by0 + cz0
        let d = plane.normal.x * plane.point.x + plane.normal.y * plane.point.y + plane.normal.z * plane.point.z;

        // Then substitue the line defined by the point to be reflected and the normal vector into the equation
        let t = 2.0 * (-d - plane.normal.x * self.x - plane.normal.y * self.y - plane.normal.z * self.z) 
            / plane.normal.magnitude_squared();
        
        let x = self.x + t * plane.normal.x;
        let y = self.y + t * plane.normal.y;
        let z = self.z + t * plane.normal.z;

        return Point{x, y, z}
    }

    pub fn reflect_through_arbitrary_plane(self, plane: Plane) -> Point {
        let displacement = self.displacement_from_plane(plane);

        return self.translate(displacement * 2.0)
    }

    /// Performs a proper rotation of this Point using a rotation Quaternion.
    /// 
    /// This method is a wrapper around the `Quaternion::rotate_point()` method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let quaternion = Quaternion::new(0.0, Vector::new(1.0, 0.0, 0.0));
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Point::new(1.0, -1.0, -1.0), point.rotate(quaternion))
    /// ```
    pub fn rotate(self, quaternion: Quaternion) -> Point {
        quaternion.rotate_point(self)
    }

    /// Performs a proper rotation Cn of this Point around a given Vector. The rotation
    /// occurs around 360/n degrees.
    /// 
    /// This method is a wrapper around the `Vector::rotate_point()` method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector = Vector::new(1.0, 0.0, 0.0);
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Point::new(1.0, -1.0, -1.0), point.rotate_around_vector(vector, 2.0))
    /// ```
    pub fn rotate_around_vector(self, vector: Vector, n: f64) -> Point {
        vector.rotate_point(self, n)
    }

    /// Performs a proper rotation Cn of this Point around a given Line. The rotation
    /// occurs around 360/n degrees.Please note that since the rotation occurs around 
    /// an axis, the position of the Line matters.
    /// 
    /// This method is a wrapper around the `Line::rotate_point()` method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let axis = Line::new(Point::new(0.0, 0.0, 1.0), Vector::new(1.0, 0.0, 0.0));
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Point::new(1.0, -1.0, 1.0), point.rotate_around_axis(axis, 2.0))
    /// ```
    pub fn rotate_around_axis(self, axis: Line, n: f64) -> Point {
        axis.rotate_point(self, n)
    }

    /// Returns a Point translated over a given Vector from the position of this Point.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point = Point::new(5.0, -1.0, 0.0);
    /// let vector = Vector::new(10.0, 2.5, -0.9);
    /// 
    /// assert_eq!(Point::new(15.0, 1.5, -0.9), point.translate(vector))
    /// ```
    pub fn translate(self, vector: Vector) -> Point {
        self + vector
    }
}



/// A vector in real 3D cartesian space, defined by its components along the x, y, and z axes.
#[derive(Debug, Copy, Clone)]
pub struct Vector {
    /// The vector component along the x-axis.
    pub x: f64,
    /// The vector component along the y-axis.
    pub y: f64,
    /// The vector component along the z-axis.
    pub z: f64
}

impl std::ops::Add<Vector> for Vector {
    type Output = Vector;

    fn add(self, rhs: Vector) -> Self::Output {
        Vector{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl std::ops::Add<Point> for Vector {
    type Output = Vector;

    fn add(self, rhs: Point) -> Self::Output {
        Vector{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl std::ops::AddAssign<Vector> for Vector {
    fn add_assign(&mut self, rhs: Vector) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl std::ops::AddAssign<Point> for Vector {
    fn add_assign(&mut self, rhs: Point) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl std::ops::Div<f64> for Vector {
    type Output = Vector;

    fn div(self, rhs: f64) -> Self::Output {
        Vector{x: self.x / rhs, y: self.y / rhs, z: self.z / rhs}
    }
}

impl ApproxEq for Vector {
    fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        return (self.x - other.x).abs() <= tolerance && (self.y - other.y).abs() <= tolerance && 
            (self.z - other.z).abs() <= tolerance
    }
}

impl std::ops::Mul<f64> for Vector {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Vector{x: self.x * rhs, y: self.y * rhs, z: self.z * rhs}
    }
}

impl std::ops::MulAssign<f64> for Vector {
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        return (self.x - other.x).abs() <= 1e-10 && (self.y - other.y).abs() <= 1e-10 && (self.z - other.z).abs() <= 1e-10
    }
}

impl VectorAlgebra<Point> for Vector {
    type Output = Vector;

    fn cross_product(self, other: Point) -> Self::Output {
        let (x, y, z) = _cross_product(self.x, self.y, self.z, other.x, other.y, other.z);

        return Vector{x, y, z}
    }
}

impl VectorAlgebra<Vector> for Vector {
    type Output = Vector;

    fn cross_product(self, other: Vector) -> Self::Output {
        let (x, y, z) = _cross_product(self.x, self.y, self.z, other.x, other.y, other.z);

        return Vector{x, y, z}
    }
}

impl Vector {
    /// Creates a new Vector from the x,y,z components.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// assert_eq!(Vector{x: 1.0, y: 1.0, z: 1.0}, Vector::new(1.0, 1.0, 1.0))
    /// ```
    pub fn new(x: f64, y: f64, z: f64) -> Vector {
        Vector{x, y, z}
    }

    /// Creates a new Vector from a Point, assuming the start point of the Vector to be [0, 0, 0].
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point = Point::new(1.5, 9.9, -5.1);
    /// assert_eq!(Vector::new(1.5, 9.9, -5.1), Vector::from_point(point))
    /// ```
    pub fn from_point(point: Point) -> Vector {
        Vector{x: point.x, y: point.y, z: point.z}
    }

    /// Creates a new Vector from a Point, assuming the start point of the Vector to be [0, 0, 0].
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point = Point::new(10.0, 0.0, 0.0);
    /// assert_eq!(Vector::new(1.0, 0.0, 0.0), Vector::from_point_normalised(point))
    /// ```
    pub fn from_point_normalised(point: Point) -> Vector {
        let magnitude = (point.x.powi(2) + point.y.powi(2) + point.z.powi(2)).sqrt();
        Vector{x: point.x / magnitude, y: point.y / magnitude, z: point.z / magnitude}
    }

    /// Creates a new Vector AB from two points A and B. The Vector starts in A and ends in B;
    /// AB = B - A.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point1 = Point::new(1.0, 0.0, -2.0);
    /// let point2 = Point::new(-8.8, -0.1, 4.3);
    /// let point3 = Point::new(2.7, 7.1, 0.9);
    /// 
    /// assert_eq!(Vector::new(-9.8, -0.1, 6.3), Vector::from_two_points(point1, point2));
    /// assert_eq!(Vector::new(9.8, 0.1, -6.3), Vector::from_two_points(point2, point1));
    /// 
    /// assert_eq!(Vector::new(1.7, 7.1, 2.9), Vector::from_two_points(point1, point3));
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), Vector::from_two_points(point1, point1))
    /// ```
    pub fn from_two_points(start: Point, end: Point) -> Vector {
        let x = end.x - start.x;
        let y = end.y - start.y;
        let z = end.z - start.z;

        return Vector{x, y, z}
    }

    /// Creates a new Vector AB from two points A and B. The Vector starts in A and ends in B;
    /// AB = B - A; but is normalised.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point1 = Point::new(1.0, 0.0, -2.0);
    /// let point2 = Point::new(-8.8, -0.1, 4.3);
    /// let point3 = Point::new(2.7, 7.1, 0.9);
    /// 
    /// assert_eq!(Vector::new(-0.841147489890757, -0.008583137651946501, 0.5407376720726295), 
    ///     Vector::from_two_points_normalised(point1, point2));
    /// assert_eq!(Vector::new(0.841147489890757, 0.008583137651946501, -0.5407376720726295), 
    ///     Vector::from_two_points_normalised(point2, point1));
    /// 
    /// assert_eq!(Vector::new(0.2164069220770561, 0.9038171451453518, 0.36916474942556626), 
    ///     Vector::from_two_points_normalised(point1, point3));
    /// 
    /// // Magnitude 0 Vectors will return a Vector of NaNs
    /// let zero_vector = Vector::from_two_points_normalised(point1, point1);
    /// assert!(zero_vector.x.is_nan());
    /// assert!(zero_vector.y.is_nan());
    /// assert!(zero_vector.z.is_nan());
    /// ```
    pub fn from_two_points_normalised(start: Point, end: Point) -> Vector {
        let x = end.x - start.x;
        let y = end.y - start.y;
        let z = end.z - start.z;

        let magnitude = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();

        return Vector{x: x / magnitude, y: y / magnitude, z: z / magnitude}
    }

    /// Returns an array containing the components of this Vector in the order of x, y, z.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector = Vector::new(1.0, -8.2, 0.0);
    /// let array = vector.as_array();
    /// 
    /// assert_eq!([1.0, -8.2, 0.0], array)
    /// ```
    pub fn as_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Returns the angle between two Vectors in radians, using the equation
    /// theta = acos(v . u / (|v| * |u|)).
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// use std::f64::consts::PI;
    /// 
    /// let x = Vector::new(1.0, 0.0, 0.0);
    /// let y = Vector::new(0.0, 1.0, 0.0);
    /// let xx = Vector::new(2.0, 0.0, 0.0);
    /// 
    /// assert_eq!(PI / 2.0, x.angle(y));
    /// assert_eq!(0.0, x.angle(xx));
    /// ```
    pub fn angle(self, other: Vector) -> f64 {
        (self.dot_product(other) / (self.magnitude() * other.magnitude())).acos()
    }

    /// Returns the angle between two Vectors in degrees, using the 
    /// Vector::angle() method and converting the result.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Vector::new(1.0, 0.0, 0.0);
    /// let y = Vector::new(0.0, 1.0, 0.0);
    /// let xx = Vector::new(2.0, 0.0, 0.0);
    /// 
    /// assert_eq!(90.0, x.angle_degrees(y));
    /// assert_eq!(0.0, x.angle_degrees(xx));
    /// ```
    pub fn angle_degrees(self, other: Vector) -> f64 {
        self.angle(other) * 180.0 / std::f64::consts::PI
    }

    /// Returns the dot product between this Vector and another Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector1 = Vector::new(1.0, 0.0, 0.0);
    /// let vector2 = Vector::new(0.0, 1.0, 0.0);
    /// let vector3 = Vector::new(0.5, -4.1, 2.0);
    /// 
    /// assert!(vector1.dot_product(vector2).abs() < 1e-10);
    /// assert!(vector2.dot_product(vector1).abs() < 1e-10);
    /// 
    /// assert!((0.5 - vector1.dot_product(vector3).abs()) < 1e-10);
    /// assert!((-4.1 - (-1.0 * vector2.dot_product(vector3)).abs()) < 1e-10);
    /// ```
    pub fn dot_product(&self, other: Vector) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Returns a Point that is the end point of this Vector in the case that a
    /// provided Point is the start point of this Vector. Contrast with the
    /// `Vector::start_point()` method which returns a start point using a 
    /// provided end point.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector = Vector::new(2.0, -1.5, 0.0);
    /// 
    /// let point1 = Point::new(0.0, 0.0, 0.0);
    /// let point2 = Point::new(1.0, 1.0, 5.0);
    /// let point3 = Point::new(-4.2, 4.2, -0.5);
    /// 
    /// assert_eq!(Point::new(2.0, -1.5, 0.0), vector.end_point(point1));
    /// assert_eq!(Point::new(3.0, -0.5, 5.0), vector.end_point(point2));
    /// assert_eq!(Point::new(-2.2, 2.7, -0.5), vector.end_point(point3));
    /// ```
    pub fn end_point(&self, start: Point) -> Point {
        Point{x: self.x + start.x, y: self.y + start.y, z: self.z + start.z}
    }

    /// Returns a Point that is the start point of this Vector in the case that a
    /// provided Point is the end point of this Vector. Contrast with the
    /// `Vector::end_point()` method which returns an end point using a provided
    /// start point.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector = Vector::new(2.0, -1.5, 0.0);
    /// 
    /// let point1 = Point::new(0.0, 0.0, 0.0);
    /// let point2 = Point::new(1.0, 1.0, 5.0);
    /// let point3 = Point::new(-4.2, 4.2, -0.5);
    /// 
    /// assert_eq!(Point::new(-2.0, 1.5, 0.0), vector.start_point(point1));
    /// assert_eq!(Point::new(-1.0, 2.5, 5.0), vector.start_point(point2));
    /// assert_eq!(Point::new(-6.2, 5.7, -0.5), vector.start_point(point3));
    /// ```
    pub fn start_point(&self, end: Point) -> Point {
        Point{x: end.x - self.x, y: end.y - self.y, z: end.z - self.z}
    }

    /// Determines whether a given Vector is a k-multiple of this Vector.
    /// This method is a wrapper around the is_approx_k_multiple()
    /// method but with the tolerance set to 1e-10.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector1 = Vector::new(1.0, 1.0, 1.0);
    /// let vector2 = Vector::new(5.0, 5.0, 5.0);
    /// let vector3 = Vector::new(1.0, 1.0, 1.5);
    /// 
    /// assert!(vector1.is_k_multiple(vector2));
    /// assert!(vector2.is_k_multiple(vector1));
    /// 
    /// assert!(!vector1.is_k_multiple(vector3));
    /// assert!(!vector2.is_k_multiple(vector3));
    /// ```
    pub fn is_k_multiple(&self, other: Vector) -> bool {
        self.is_approx_k_multiple(other, 1e-10)
    }

    /// Determines whether a given Vector is a k-multiple of this Vector, 
    /// taking into account the provided tolerance. The two Vectors are
    /// k multiples if the qutient of each set of coordinates of the Vectors
    /// is equal within the tolerance.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector1 = Vector::new(1.0, 1.0, 1.0);
    /// let vector2 = Vector::new(5.0, 5.0, 5.0);
    /// let vector3 = Vector::new(1.0, 1.0, 1.5);
    /// 
    /// assert!(vector1.is_approx_k_multiple(vector2, 1e-10));
    /// assert!(vector2.is_approx_k_multiple(vector1, 1e-1));
    /// 
    /// assert!(!vector1.is_approx_k_multiple(vector3, 1e-10));
    /// assert!(!vector2.is_approx_k_multiple(vector3, 1e-1));
    /// 
    /// assert!(vector1.is_approx_k_multiple(vector3, 10.0));
    /// ```
    pub fn is_approx_k_multiple(&self, other: Vector, tolerance: f64) -> bool {
        // First check if there is no mismatch between zero values, and if not, compute the ratio of the vectors
        let mut result = [0.0; 3];
        for (i, (val1, val2)) in core::iter::zip(self.as_array(), other.as_array()).enumerate() {
            if (val1.abs() < 1e-10 && val2.abs() > 1e-10) || (val1.abs() > 1e-10 && val2.abs() < 1e-10) {
                if (val1 - val2).abs() < tolerance {
                    result[i] = 0.0 / 0.0
                } else {
                    return false
                }
            } else if val1.abs() < tolerance && val2.abs() < tolerance {
                result[i] = 0.0 / 0.0
            } else {
                result[i] = val1 / val2
            }
        };

        // Then, check if the ratio is the same for all the parameters
        if (result[0] - result[1]).abs() > tolerance || (result[0] - result[2]).abs() > tolerance || (result[1] - result[2]).abs() > tolerance {
            return false
        } else {
            return true
        }
    }

    /// Returns the magnitude (i.e. the length) of this Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector1 = Vector::new(1.0, 0.0, 0.0);
    /// let vector2 = Vector::new(6.0, 2.0, -3.0);
    /// 
    /// assert!((1.0 - vector1.magnitude()).abs() <= 1e-10);
    /// assert!((7.0 - vector2.magnitude()).abs() <= 1e-10)
    /// ```
    pub fn magnitude(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    /// Returns the square of the magnitude (i.e. the length) of this Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector1 = Vector::new(1.0, 0.0, 0.0);
    /// let vector2 = Vector::new(6.0, 2.0, -3.0);
    /// 
    /// assert!((1.0 - vector1.magnitude_squared()).abs() <= 1e-10);
    /// assert!((49.0 - vector2.magnitude_squared()).abs() <= 1e-10)
    /// ```
    pub fn magnitude_squared(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }

    pub fn normalise(self) -> Vector {
        self / self.magnitude()
    }

    /// Returns a Vector that is in the same position as if this Vector were
    /// rotated over 360/n degrees around the provided Vector.
    /// 
    /// In other words, the Vector that this method is called on is rotated
    /// around the Vector that is passed to the method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector_to_rotate = Vector::new(1.0, 0.0, 0.0);
    /// let rotation_vector = Vector::new(1.0, 1.0, 0.0);
    /// 
    /// assert_eq!(Vector::new(0.0, 1.0, 0.0), vector_to_rotate.rotate_around(rotation_vector, 2.0))
    /// ```
    pub fn rotate_around(self, other: Vector, n: f64) -> Vector {
        let quaternion = Quaternion::create_rotation_quaternion(other, n);

        return quaternion.rotate_vector(self)
    }

    /// Returns a Point that is in the same position as if the given Point were
    /// rotated over 360/n degrees around this Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let vector = Vector::new(1.0, 0.0, 0.0);
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Point::new(1.0, -1.0, -1.0), vector.rotate_point(point, 2.0))
    /// ```
    pub fn rotate_point(self, point: Point, n: f64) -> Point {
        let quaternion = Quaternion::create_rotation_quaternion(self, n);

        return quaternion.rotate_point(point)
    }

    /// Returns a Point that is in the same position as if the given Point were
    /// rotated over 360/n degrees around this Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let rotation_vector = Vector::new(1.0, 0.0, 0.0);
    /// let vector_to_rotate = Vector::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Vector::new(1.0, -1.0, -1.0), rotation_vector.rotate_vector(vector_to_rotate, 2.0))
    /// ```
    pub fn rotate_vector(self, other: Vector, n: f64) -> Vector {
        let quaternion = Quaternion::create_rotation_quaternion(self, n);

        return quaternion.rotate_vector(other)
    }
}



/// A quaternion q + xi + yj + zk, where i,j,k are basis vectors.
#[derive(Debug, Copy, Clone)]
pub struct Quaternion {
    /// The constant coefficient, q0.
    pub q: f64,
    /// A Vector representing the quaternion coefficients that correspond to basis vectors, q1*i + q2*j + q3.k
    pub vector: Vector
}

impl ApproxEq for Quaternion {
    fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        return (self.q - other.q).abs() <= tolerance && self.vector.approx_eq(&other.vector, tolerance)
    }
}

impl PartialEq for Quaternion {
    fn eq(&self, other: &Self) -> bool {
        (self.q - other.q).abs() <= 1e-10 && self.vector == other.vector
    }
}

impl Quaternion {
    /// Creates a new Quaternion from the constant coefficient and a Vector that corresponds to vector coefficients.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let quaternion = Quaternion::new(0.0, Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert_eq!(Quaternion{q: 0.0, vector: Vector::new(1.0, 0.0, 0.0)}, quaternion)
    /// ```
    pub fn new(q: f64, vector: Vector) -> Quaternion {
        Quaternion {q, vector}
    }

    /// Creates a new Quaternion such that it can represent a rotation around the provided Vector around
    /// 360/n degrees. The following equation is used for this, where theta is half of the angle being
    /// rotated through and vector is the vector being rotated around:
    /// q = cos(theta/2) + (vector/|vector|)*sin(theta/2)
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let quaternion = Quaternion::create_rotation_quaternion(Vector::new(1.0, 0.0, 0.0), 2.0);
    /// 
    /// assert_eq!(Quaternion::new(0.0, Vector::new(1.0, 0.0, 0.0)), quaternion)
    /// ```
    pub fn create_rotation_quaternion(vector: Vector, n: f64) -> Quaternion {
        // Angle = 2/n * pi → theta/2 = 1/n*pi

        let angle = 1.0 / n  * std::f64::consts::PI;
        let magnitude = vector.magnitude();
        let sine = angle.sin();

        let q = angle.cos();
        let x = vector.x / magnitude * sine;
        let y = vector.y / magnitude * sine;
        let z = vector.z / magnitude * sine;

        return Quaternion{q, vector: Vector{x, y, z}}
    }

    /// Returns a Point that is in such a position as if it were rotated using this
    /// Quaternion. In other words, rotates the provided Point using this Quaternion and returns a
    /// new Point.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let quaternion = Quaternion::create_rotation_quaternion(Vector::new(1.0, 0.0, 0.0), 2.0);
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Point::new(1.0, -1.0, -1.0), quaternion.rotate_point(point))
    /// ```
    pub fn rotate_point(&self, point: Point) -> Point {
        // Compute parameters for rotation using qvq* = av + bq + c(q×v), where a = q0^2 - |q|^2, b = 2(v.q), c = 2q0
        let a: f64 = self.q.powi(2) - (self.vector.x.powi(2) + self.vector.y.powi(2) + self.vector.z.powi(2));
        let b = 2.0 * (point.x * self.vector.x + point.y * self.vector.y + point.z * self.vector.z);
        let c = 2.0 * self.q;
        let cross = self.vector.cross_product(point);

        let x = a * point.x + b * self.vector.x + c * cross.x;
        let y = a * point.y + b * self.vector.y + c * cross.y;
        let z = a * point.z + b * self.vector.z + c * cross.z;

        return Point {x, y, z}
    }

    /// Returns a Vector that is in such a position as if it were rotated using this
    /// Quaternion. In other words, rotates the provided Vector using this Quaternion and returns a
    /// new Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let quaternion = Quaternion::create_rotation_quaternion(Vector::new(1.0, 0.0, 0.0), 2.0);
    /// let vector = Vector::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Vector::new(1.0, -1.0, -1.0), quaternion.rotate_vector(vector))
    /// ```
    pub fn rotate_vector(&self, vector: Vector) -> Vector {
        // Compute parameters for rotation using qvq* = av + bq + c(q×v), where a = q0^2 - |q|^2, b = 2(v.q), c = 2q0
        let a: f64 = self.q.powi(2) - (self.vector.x.powi(2) + self.vector.y.powi(2) + self.vector.z.powi(2));
        let b = 2.0 * (vector.x * self.vector.x + vector.y * self.vector.y + vector.z * self.vector.z);
        let c = 2.0 * self.q;
        let cross = self.vector.cross_product(vector);

        let x = a * vector.x + b * self.vector.x + c * cross.x;
        let y = a * vector.y + b * self.vector.y + c * cross.y;
        let z = a * vector.z + b * self.vector.z + c * cross.z;

        return Vector {x, y, z}
    }
}



/// A physical Line in real 3D space, defined by a Point that lies on the Line and a Vector
/// signifying the direction of the Line.
#[derive(Debug, Copy, Clone)]
pub struct Line{
    /// A Point that lies on the Line
    pub point: Point,
    /// A Vector which determines the direction of the Line.
    pub vector: Vector
}

impl ApproxEq for Line {
    fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        return self.vector.is_approx_k_multiple(other.vector, tolerance) && other.point.is_approx_on_line(*self, tolerance)
    }
}

impl PartialEq for Line {
    fn eq(&self, other: &Self) -> bool {
        return self.vector.is_k_multiple(other.vector) && other.point.is_on_line(*self)
    }
}

impl Line {
    /// Creates a new Line from a Point and a Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// let result = Line{point: Point::new(0.0, 0.0, 0.0), vector: Vector::new(1.0, 0.0, 0.0)};
    /// 
    /// assert_eq!(result, line)
    /// ```
    pub fn new(point: Point, vector: Vector) -> Line {
        Line{point, vector}
    }

    /// Creates a new Line using two Points to define it.
    /// A Line is created only if the provided Points are not
    /// identical, otherwise None is returned.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let valid_line = Line::from_points(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0));
    /// assert_eq!(Some(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0))), valid_line);
    /// 
    /// let invalid_line = Line::from_points(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0));
    /// assert_eq!(None, invalid_line);
    /// ```
    pub fn from_points(start: Point, end: Point) -> Option<Line> {
        if start == end {
            return None
        } else {
            Some(Line{point: start, vector: Vector::from_two_points(start, end)})
        }
    }

    /// Creates a new Line using two Points to define it,
    /// with the directional Vector normalised.
    /// A Line is created only if the provided Points are not
    /// identical, otherwise None is returned.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let valid_line = Line::from_points_normalised(Point::new(0.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0));
    /// assert_eq!(Some(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0))), valid_line);
    /// 
    /// let invalid_line = Line::from_points_normalised(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0));
    /// assert_eq!(None, invalid_line);
    /// ```
    pub fn from_points_normalised(start: Point, end: Point) -> Option<Line> {
        if start == end {
            return None
        } else {
            Some(Line{point: start, vector: Vector::from_two_points_normalised(start, end)})
        }
    }

    /// Creates a new Line using two Points to define it.
    /// A Line is created only if the provided Points are
    /// approximately different (determined using Point::approx_eq()), 
    /// otherwise None is returned.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let valid_line = Line::from_points_tolerance(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), 1e-5);
    /// assert_eq!(Some(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0))), valid_line);
    /// 
    /// let invalid_line = Line::from_points_tolerance(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0), 1e-5);
    /// assert_eq!(None, invalid_line);
    /// 
    /// let points_within_tolerance = Line::from_points_tolerance(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), 10.0);
    /// assert_eq!(None, points_within_tolerance);
    /// ```
    pub fn from_points_tolerance(start: Point, end: Point, tolerance: f64) -> Option<Line> {
        if start.approx_eq(&end, tolerance) {
            return None
        } else {
            Some(Line{point: start, vector: Vector::from_two_points(start, end)})
        }
    }

    /// Creates a new Line using two Points to define it.
    /// The directional Vector of the created line will be
    /// normalised.
    /// 
    /// A Line is created only if the provided Points are
    /// approximately different (determined using Point::approx_eq()), 
    /// otherwise None is returned.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let valid_line = Line::from_points_tolerance_normalised(Point::new(0.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0), 1e-5);
    /// assert_eq!(Some(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0))), valid_line);
    /// 
    /// let invalid_line = Line::from_points_tolerance_normalised(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 0.0), 1e-5);
    /// assert_eq!(None, invalid_line);
    /// 
    /// let points_within_tolerance = Line::from_points_tolerance_normalised(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), 10.0);
    /// assert_eq!(None, points_within_tolerance);
    /// ```
    pub fn from_points_tolerance_normalised(start: Point, end: Point, tolerance: f64) -> Option<Line> {
        if start.approx_eq(&end, tolerance) {
            return None
        } else {
            Some(Line{point: start, vector: Vector::from_two_points_normalised(start, end)})
        }
    }

    /// Computes the angle in radians between two Lines, regardless
    /// of their mutual position. This is done by computing
    /// the angle between the directional Vectors of the
    /// two lines.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// use std::f64::consts::PI;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// let y = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0));
    /// let xx = Line::new(Point::new(0.0, 0.0, 10.0), Vector::new(2.0, 0.0, 0.0));
    /// 
    /// assert_eq!(PI / 2.0, x.angle(y));
    /// assert_eq!(0.0, x.angle(xx));
    /// ```
    pub fn angle(self, other: Line) -> f64 {
        self.vector.angle(other.vector)
    }

    /// Computes the angle in degrees between two Lines, regardless
    /// of their mutual position. This is done by computing
    /// the angle between the directional Vectors of the
    /// two lines.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// use std::f64::consts::PI;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// let y = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0));
    /// let xx = Line::new(Point::new(0.0, 0.0, 10.0), Vector::new(2.0, 0.0, 0.0));
    /// 
    /// assert_eq!(90.0, x.angle_degrees(y));
    /// assert_eq!(0.0, x.angle_degrees(xx));
    /// ```
    pub fn angle_degrees(self, other: Line) -> f64 {
        self.vector.angle_degrees(other.vector)
    }

    pub fn fast_approx_eq(self, other: Line, tolerance: f64) -> bool {
        return self.vector.is_approx_k_multiple(other.vector, tolerance)
    }

    /// Determines whether a given Point lies on this Line. This is a wrapper around the
    /// `Point::lies_on_line()` method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert!(line.has_point(Point::new(0.0, 0.0, 0.0)));
    /// assert!(line.has_point(Point::new(5.9, 0.0, 0.0)));
    /// assert!(line.has_point(Point::new(-0.2, 0.0, 0.0)));
    /// assert!(!line.has_point(Point::new(0.0, 0.5, 0.0)))
    /// ```
    pub fn has_point(self, point: Point) -> bool {
        point.is_on_line(self)
    }

    /// Determines whether a given Point lies on this Line. This is a wrapper around the
    /// `Point::lies_on_line()` method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert!(Point::new(0.0, 0.0, 0.0).is_approx_on_line(line, 1e-10));
    /// assert!(Point::new(5.9, 0.0, 0.0).is_approx_on_line(line, 1e-5));
    /// assert!(Point::new(-0.2, 0.0, 0.0).is_approx_on_line(line, 1e-1));
    /// assert!(!Point::new(0.0, 0.5, 0.0).is_approx_on_line(line, 1e-10));
    /// 
    /// assert!(Point::new(1.0, 0.00001, 0.0).is_approx_on_line(line, 1e-4));
    /// assert!(!Point::new(1.0, 0.00001, 0.0).is_approx_on_line(line, 1e-6));
    /// ```
    pub fn approx_has_point(self, point: Point, tolerance: f64) -> bool {
        point.is_approx_on_line(self, tolerance)
    }

    /// Computes a Point that lies on this line, given a factor t. It does this by substituting t into the
    /// parametric equations of a line.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert_eq!(Point::new(1.0, 0.0, 0.0), line.compute_point(1.0));
    /// assert_eq!(Point::new(5.0, 0.0, 0.0), line.compute_point(5.0));
    /// assert_eq!(Point::new(-3.14, 0.0, 0.0), line.compute_point(-3.14));
    /// ```
    pub fn compute_point(self, t: f64) -> Point {
        let x = self.point.x + self.vector.x * t;
        let y = self.point.y + self.vector.y * t;
        let z = self.point.z + self.vector.z * t;

        Point{x, y, z}
    }

    pub fn normalise(self) -> Line {
        Line{point: self.point, vector: self.vector.normalise()}
    }

    /// Returns a Point such that its position is the same as if the given Point were rotated
    /// 360/n degrees around this Line.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let axis = Line::new(Point::new(0.0, 0.0, 1.0), Vector::new(1.0, 0.0, 0.0));
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// 
    /// assert_eq!(Point::new(1.0, -1.0, 1.0), point.rotate_around_axis(axis, 2.0))
    /// ```
    pub fn rotate_point(self, point: Point, n: f64) -> Point {
        let t = -1.0 * self.point.dot_product(self.vector) / self.vector.magnitude_squared();
        
        let displacement = Vector::from_two_points(self.compute_point(t), Point{x: 0.0, y: 0.0, z: 0.0});

        return point.translate(displacement).rotate_around_vector(self.vector, n).translate(displacement * -1.0)
    }

    /// Returns a Vec of Points where each Point is rotated 360/n degrees around this Line.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let axis = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// let coordinates = vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0), Point::new(0.0, 0.0, 1.0)];
    /// 
    /// assert_eq!(vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, -1.0, 0.0), Point::new(0.0, 0.0, -1.0)], 
    ///     axis.rotate_coordinates(coordinates, 2.0, 3))
    /// ```
    pub fn rotate_coordinates<T: IntoIterator<Item = Point>>(self, coordinates: T, n: f64, length: usize) -> Vec<Point> {
        let mut result = Vec::with_capacity(length);

        for point in coordinates {
            result.push(self.rotate_point(point, n))
        }

        return result
    }
}



/// A Plane in real 3D space as defined by a Point that lies on the Plane and
/// a Vector that is normal to the Plane.
#[derive(Debug, Copy, Clone)]
pub struct Plane {
    /// A Point that lies on the Plane.
    pub point: Point,
    /// The Vector normal to the Plane.
    pub normal: Vector
}

impl ApproxEq for Plane {
    fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        self.approx_has_point(other.point, tolerance) && self.normal.is_approx_k_multiple(other.normal, tolerance)
    }
}

impl PartialEq for Plane {
    fn eq(&self, other: &Self) -> bool {
        self.has_point(other.point) && self.normal.is_k_multiple(other.normal)
    }
}

impl Plane {
    /// Creates a new Plane from a Point and a Vector so long as the Vector is non-zero.
    /// If the Vector is quivalent to Vector{x: 0.0, y: 0.0, z: 0.0}, None is returned
    /// since a Plane must be defined by a valid Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0));
    /// let result = Plane{point: Point::new(0.0, 0.0, 0.0), normal: Vector::new(0.0, 1.0, 0.0)};
    /// 
    /// assert_eq!(Some(result), plane);
    /// assert_eq!(None, Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, -0.0, 0.0)));
    /// ```
    pub fn new(point: Point, normal: Vector) -> Option<Plane> {
        if normal.x.abs() < 1e-10 && normal.y.abs() < 1e-10 && normal.z.abs() < 1e-10 {
            return None
        } else {
            return Some(Plane{point, normal})
        }
    }

    /// Creates a new Plane from three points if it is possible. If the two vectors resulting
    /// from the three points are parallel, no Plane can be created and so None is returned.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point1 = Point::new(0.0, 0.0, 0.0);
    /// let point2 = Point::new(1.0, 0.0, 0.0);
    /// let point3 = Point::new(0.0, 0.0, 1.0);
    /// 
    /// let plane = Plane::from_three_points(point1, point2, point3);
    /// 
    /// assert_eq!(Some(Plane::new(point1, Vector::new(0.0, -1.0, 0.0)).unwrap()), plane);
    /// assert_eq!(None, Plane::from_three_points(point1, point2, Point::new(-1.0, 0.0, 0.0)))
    /// ```
    pub fn from_three_points(point1: Point, point2: Point, point3: Point) -> Option<Plane> {
        let vector1 = Vector::from_two_points(point1, point2);
        let vector2 = Vector::from_two_points(point1, point3);

        let normal = vector1.cross_product(vector2);

        if normal.x.abs() < 1e-10 && normal.y.abs() < 1e-10 && normal.z.abs() < 1e-10 {
            return None
        } else {
            return Some(Plane{point: point1, normal})
        }
    }

    /// Creates a new Plane from three points if it is possible. If the two vectors resulting
    /// from the three points are parallel, no Plane can be created and so None is returned.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point1 = Point::new(0.0, 0.0, 0.0);
    /// let point2 = Point::new(1.0, 0.0, 0.0);
    /// let point3 = Point::new(0.0, 0.0, 1.0);
    /// 
    /// let plane = Plane::from_three_points_tolerance(point1, point2, point3, 1e-10);
    /// 
    /// assert_eq!(Some(Plane::new(point1, Vector::new(0.0, -1.0, 0.0)).unwrap()), plane);
    /// assert_eq!(None, Plane::from_three_points_tolerance(point1, point2, Point::new(-1.0, 0.0, 0.0), 1e-10))
    /// ```
    pub fn from_three_points_tolerance(point1: Point, point2: Point, point3: Point, tolerance: f64) -> Option<Plane> {
        let vector1 = Vector::from_two_points(point1, point2);
        let vector2 = Vector::from_two_points(point1, point3);

        let normal = vector1.cross_product(vector2);

        if normal.x.abs() < tolerance && normal.y.abs() < tolerance && normal.z.abs() < tolerance {
            return None
        } else {
            return Some(Plane{point: point1, normal})
        }
    }

    pub fn from_multiple_points(points: &Vec<&Point>, tolerance: f64) -> Option<Plane> {        
        let length = points.len() as f64;

        if length < 3.0 {
            return None
        }

        let mut centroid = Point::new(0.0, 0.0, 0.0);
        for point in points {
            centroid += *point
        }
        centroid = centroid / length;

        // Calculate full 3x3 covariance matrix, excluding symmetries:
        let mut xx = 0.0; let mut xy = 0.0; let mut xz = 0.0;
        let mut yy = 0.0; let mut yz = 0.0; let mut zz = 0.0;

        for p in points {
            let r = *p - centroid;
            xx += r.x * r.x;
            xy += r.x * r.y;
            xz += r.x * r.z;
            yy += r.y * r.y;
            yz += r.y * r.z;
            zz += r.z * r.z;
        }

        xx /= length; xy /= length; xz /= length;
        yy /= length; yz /= length; zz /= length;

        let mut weighted_dir = Vector::new(0.0, 0.0, 0.0);
        for axis in [0, 1, 2] {
            weighted_dir += _linear_regression(axis, xx, xy, xz, yy, yz, zz, weighted_dir)
                .expect("Axis should never reach a value different from 0, 1, 2")
        }

        weighted_dir = weighted_dir.normalise();

        if weighted_dir.x.abs() > tolerance || weighted_dir.y.abs() > tolerance || weighted_dir.z.abs() > tolerance {
            return Some(Plane{point: centroid, normal: weighted_dir})
        } else {
            return None
        }
    }

    pub fn from_multiple_points2(points: &Vec<&Point>, tolerance: f64) -> Option<Plane> {
        let length = points.len();
        let mut vectors = Vec::with_capacity(length * (length-1));

        let mut centre = Point::new(0.0, 0.0, 0.0);
        for (i, p1) in points.iter().enumerate() {
            centre += **p1;
            for p2 in points[i+1..].iter() {
                vectors.push(Vector::from_two_points_normalised(**p1, **p2));
            }
        }
        centre = centre / length as f64;

        let mut normal = Vector::new(0.0, 0.0, 0.0);
        for (i, v1) in vectors.iter().enumerate() {
            for v2 in vectors[i+1..].iter() {
                normal += v1.cross_product(*v2)
            }
        }
        normal = normal / length as f64;

        if normal.x.abs() > tolerance || normal.y.abs() > tolerance || normal.z.abs() > tolerance {
            return Some(Plane{point: centre, normal: normal.normalise()})
        } else {
            return None
        }
    }

    /// Returns a Point such that its position is the same as if it were reflected
    /// through this Plane.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let point = Point::new(1.0, 1.0, 1.0);
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(Point::new(1.0, -1.0, 1.0), point.reflect(plane))
    /// ```
    pub fn reflect_point(self, point: Point) -> Point {
        point.reflect_through_arbitrary_plane(self)
    }

    /// Determines whether a given Point lies on this Plane.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
    /// 
    /// assert!(Point::new(0.0, 0.0, 0.0).is_on_plane(plane));
    /// assert!(Point::new(9.5, 0.0, 0.0).is_on_plane(plane));
    /// assert!(Point::new(0.0, 0.0, -1.8).is_on_plane(plane));
    /// assert!(Point::new(-11.0, 0.0, 5.6).is_on_plane(plane));
    /// 
    /// assert!(!Point::new(0.0, 1.4, 0.0).is_on_plane(plane))
    /// ```
    pub fn has_point(self, point: Point) -> bool {
        point.is_on_plane(self)
    }

    /// Determines whether a given Point lies on this Plane.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)).unwrap();
    /// 
    /// assert!(Point::new(0.0, 0.0, 0.0).is_on_plane(plane));
    /// assert!(Point::new(9.5, 0.0, 0.0).is_on_plane(plane));
    /// assert!(Point::new(0.0, 0.0, -1.8).is_on_plane(plane));
    /// assert!(Point::new(-11.0, 0.0, 5.6).is_on_plane(plane));
    /// 
    /// assert!(!Point::new(0.0, 1.4, 0.0).is_on_plane(plane))
    /// ```
    pub fn approx_has_point(self, point: Point, tolerance: f64) -> bool {
        point.is_approx_on_plane(self, tolerance)
    }

    /// Computes the displacement of this Plane from a given Point.
    /// I.e., it computes the Vector through which this Plane would
    /// have to be translated so that the given Point lies on this
    /// Plane.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, 2.0)).unwrap();
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), plane.displacement_from_point(Point::new(0.0, 0.0, 5.0)));
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), plane.displacement_from_point(Point::new(5.0, 5.0, 5.0)));
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, 1.0), plane.displacement_from_point(Point::new(0.0, 0.0, 6.0)));
    /// assert_eq!(Vector::new(0.0, 0.0, -1.0), plane.displacement_from_point(Point::new(0.0, 0.0, 4.0)));
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, 1.0), plane.displacement_from_point(Point::new(5.0, 0.0, 6.0)));
    /// assert_eq!(Vector::new(0.0, 0.0, -1.0), plane.displacement_from_point(Point::new(0.0, 5.0, 4.0)));
    /// ```
    pub fn displacement_from_point(&self, point: Point) -> Vector {
        let distance = point.distance_from_plane(*self);
        return self.normal * distance / self.normal.magnitude()
    }

    /// Computes the d parameter in the ax + by + cz + d = 0 equation
    /// of a Plane. In this equation, (a, b, c) is the normal Vector,
    /// (x, y, z) is any Point on this Plane, and d is the negative
    /// dot product between the normal Vector and the Point defining
    /// this plane.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// assert_eq!(-5.0, Plane::new(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, 1.0)).unwrap().d());
    /// assert_eq!(0.0, Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap().d());
    /// ```
    pub fn d(&self) -> f64 {
        -1.0 * (self.normal.x * self.point.x + self.normal.y * self.point.y + self.normal.z * self.point.z)
    }
}


fn _linear_regression(axis: u16, xx: f64, xy: f64, xz: f64, yy: f64, yz: f64, zz: f64,
    weighted_dir: Vector) -> Result<Vector, ()> {
    let (det, x, y, z) = if axis == 0 {
        let det = yy*zz - yz*yz;
        (det, det, xz*yz - xy*zz, xy*yz - xz*yy)
    } else if axis == 1 {
        let det = xx*zz - xz*xz;
        (det, xz*yz - xy*zz, det, xy*xz - yz*xx)
    } else if axis == 2 {
        let det = xx*yy - xy*xy;
        (det, xy*yz - xz*yy, xy*xz - yz*xx, det)
    } else {
        return Err(())
    };

    let axis_dir = Vector::new(x, y, z);

    let weight = if weighted_dir.dot_product(axis_dir) < 0.0 {
        - det * det
    } else {
        det * det
    };

    return Ok(axis_dir * weight)
}





#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[test]
    fn test_reflect_point_through_arbitrary_plane() {
        let plane = Plane::new(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, 1.0)).unwrap();

        assert_eq!(Point::new(0.0, 0.0, 4.0), Point::new(0.0, 0.0, 6.0).reflect_through_arbitrary_plane(plane));
        assert_eq!(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 10.0).reflect_through_arbitrary_plane(plane));
        assert_eq!(Point::new(0.0, 0.0, 15.0), Point::new(0.0, 0.0, -5.0).reflect_through_arbitrary_plane(plane));

        let plane = Plane::new(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, -1.0)).unwrap();

        assert_eq!(Point::new(0.0, 0.0, 4.0), Point::new(0.0, 0.0, 6.0).reflect_through_arbitrary_plane(plane));
        assert_eq!(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 10.0).reflect_through_arbitrary_plane(plane));
        assert_eq!(Point::new(0.0, 0.0, 15.0), Point::new(0.0, 0.0, -5.0).reflect_through_arbitrary_plane(plane));
    }

    #[test]
    fn test_reflect_point_origin() {
        let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();

        assert_eq!(Point::new(0.0, 0.0, -1.0), Point::new(0.0, 0.0, 1.0).reflect_through_arbitrary_plane(plane))
    }

    #[test]
    fn test_rotate_point_around_vector_x_axis_180() {
        let point = Point::new(1.0, 1.0, 1.0);
        let axis = Vector::new(1.0, 0.0, 0.0);
        let result = axis.rotate_point(point, 2.0);

        assert_eq!(Point::new(1.0, -1.0, -1.0), result)
    }

    #[test]
    fn test_rotate_point_around_vector_point_on_axis_180() {
        let point = Point::new(0.5, 0.0, 0.0);
        let axis = Vector::new(1.0, 0.0, 0.0);
        let result = axis.rotate_point(point, 2.0);

        assert_eq!(result, point)
    }
    // TODO: check tests with external rotation source
    #[test]
    fn test_rotate_point_around_vector_3d_axis_180() {
        let point = Point::new(-1.0, 2.0, 5.0);
        let axis = Vector::new(1.0, 1.0, 1.0);
        let result = axis.rotate_point(point, 2.0);

        assert_eq!(Point::new(5.0, 2.0, -1.0), result)
    }

    #[test]
    fn test_rotate_point_around_vector_x_axis_90() {
        let point = Point::new(1.0, 1.0, 0.0);
        let axis = Vector::new(1.0, 0.0, 0.0);
        let result = axis.rotate_point(point, 4.0);

        assert_eq!(Point::new(1.0, 0.0, 1.0), result)
    }

    #[test]
    fn test_rotate_point_around_vector_x_axis_120() {
        let point = Point::new(1.0, 1.0, 0.0);
        let axis = Vector::new(1.0, 0.0, 0.0);
        let result = axis.rotate_point(point, 3.0);

        assert_eq!(Point::new(1.0, -0.5, 0.8660254037844388), result)
    }

    #[rstest]
    #[case::x_180(Point::new(1.0, 1.0, 1.0), Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)), 
        2.0, Point::new(1.0, -1.0, -1.0))]
    #[case::on_line_180(Point::new(0.5, 0.0, 0.0), Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        2.0, Point::new(0.5, 0.0, 0.0))]
    #[case::line_3d_180(Point::new(-1.0, 2.0, 5.0), Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 1.0, 1.0)),
        2.0, Point::new(5.0, 2.0, -1.0))]
    #[case::x_90(Point::new(1.0, 1.0, 0.0), Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        4.0, Point::new(1.0, 0.0, 1.0))]
    #[case::x_120(Point::new(1.0, 1.0, 0.0), Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        3.0, Point::new(1.0, -0.5, 0.8660254037844388))]
    #[case::parallel_x_180(Point::new(1.0, 1.0, 0.0), Line::new(Point::new(0.0, 2.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        2.0, Point::new(1.0, 3.0, 0.0))]
    #[case::parallel_y_180(Point::new(-0.05, 11.41, 2.27), Line::new(Point::new(-0.05, 10.15, 2.27),
        Vector::new(1e-10, 1.25, 1e-10)), 2.0, Point::new(-0.04999999979840001, 11.41, 2.2700000002016005))]
    fn test_rotate_point_around_axis(#[case] point: Point, #[case] line: Line, #[case] n: f64, #[case] expected: Point) {
        let result = line.rotate_point(point, n);
        assert_eq!(expected, result)
    }
}