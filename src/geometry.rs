//! 


/// Used for determining approximate equality between two objects.
/// 
/// Therefore, allows for two objects to be equal if they are within
/// certain tolerance. Rust's type system ensures that the two
/// objects being compared are of the same type, but the details
/// of that equality is determined by the implementation of this trait.
/// 
/// # Examples
/// 
/// ```
/// use point_group::geometry::ApproxEq;
/// 
/// #[derive(Copy, Clone)]
/// struct Point {
///     x: f64,
///     y: f64,
/// }
/// 
/// impl ApproxEq for Point {
///     fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
///         (self.x - other.x).abs() < tolerance && (self.y - other.y).abs() < tolerance
///     }
/// }
/// 
/// let point = Point{x: 0.0, y: 0.0};
/// 
/// assert!(point.approx_eq(&Point{x: 0.0, y: 0.0}, 1e-10));
/// assert!(!point.approx_eq(&Point{x: 1.0, y: 0.0}, 1e-10));
/// assert!(point.approx_eq(&Point{x: 1.0, y: 0.0}, 10.0));
/// ```
pub trait ApproxEq {
    fn approx_eq(&self, other: &Self, tolerance: f64) -> bool;
}



/// Used for implementing vector algebra on vector-like objects,
/// in this case the [`Point`] struct and [`Vector`] struct.
/// 
/// In reality, this trait is just an abstraction used for 
/// overloading the cross_product and dot_product methods
/// such that they can accept and output both [`Point`] and
/// [`Vector`]. The actual math implementation is handled by
/// private functions.
pub trait VectorAlgebra<T> {
    type Output;

    fn cross_product(self, other: T) -> Self::Output;

    fn dot_product(self, other: T) -> f64;
}

fn _cross_product(rx: f64, ry: f64, rz: f64, lx: f64, ly: f64, lz: f64) -> (f64, f64, f64) {
    let x = ry * lz - rz * ly;
    let y = rz * lx - rx * lz;
    let z = rx * ly - ry * lx;

    return (x, y, z)
}

fn _dot_product(rx: f64, ry: f64, rz: f64, lx: f64, ly: f64, lz: f64) -> f64 {
    rx * lx + ry * ly + rz * lz
}



/// Used for determining the relative positions
/// of two geometric objects.
pub trait RelativePositionTrait<T> {
    /// Used for computing the displacement between two
    /// geometric objects. The displacement is always
    /// expressed in the form of a [`Vector`] which
    /// points in the direction of the object on which
    /// this method is called to the the other object.
    fn displacement(self, other: T) -> Vector;

    /// Used for computing the distance between two 
    /// geometric objects.
    fn distance(self, other: T) -> f64;

    /// Used for determining the relative positions of two
    /// geometric objects, which is expressed in the form
    /// of the variants of the [`RelativePosition`] enum.
    /// 
    /// The main use case of this method for when various 
    /// properties, such as displacement or distance, are
    /// calculated in a different way depending on the 
    /// relative position of the two objects.
    fn relative_position(self, other: T) -> RelativePosition;
}


/// Used for modelling the relative position of two geometric
/// objects.
/// 
/// **Note:** Some combinations of objects can adopt only some
/// but not all of these relative positions, or, more accurately,
/// some of the variants represent the same configuration. For 
/// example, two [`Point`]s can only have two configurations: 
/// either they are exactly equal or they are not. In these cases,
/// the priority in nomenclature is Equivalent before Intersecting 
/// and Parallel before Skew. As such, the two configurations of
/// two Points are called `Equivalent` and `Parallel` respectively.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum RelativePosition {
    /// Represents the equivalence of two objects. In regards to 
    /// this enum, two objects are `Equivalent` when they occupy the
    /// same space, which can take two forms:
    ///  - the two objects are exactly equal
    ///  - one of the objects is a part of the other objects, or
    /// more mathematically, the set of Points of one of the objects
    /// is a subset of the set of Points of the other object.
    /// 
    /// Contrast this with the `Intersecting` variant, where the
    /// two objects share exactly one Point. **Note:** in cases
    /// where Equivalent and Intersecting are both equally valid
    /// descriptors, such as for two [`Point`]s, the `Equivalent`
    /// variant takes higher priority, and so only `Equivalent` is
    /// used to describe the configuration.
    Equivalent,

    /// Represents the configuration of two objects where they
    /// are parallel to one another. In regards to this enum, 
    /// two objects are considered `Parallel` when they do not
    /// share any [`Point`]s and their directions are the same.
    /// 
    /// Contrast this with the `Skew` variant, where the objects 
    /// also don't share any Points but where their directions 
    /// are different. **Note:** however, where the relative
    /// direction of the two objects is not relevant, either
    /// because they have no directions (like two Points) or
    /// because they can only have one relative direction when
    /// not sharing any Points (like a [`Line`] and a [`Plane`]),
    /// the `Parallel` variant takes priority over `Skew`, and
    /// is the only one used to describe such configurations.
    Parallel,

    /// Reperesents such positioning of two objects where they 
    /// share exactly one [`Point`].
    /// 
    /// Contrast this with the `Equivalent` variant, where they
    /// share the same space. **Note:** in cases where Equivalent 
    /// and Intersecting are both equally valid descriptors, such 
    /// as for two [`Point`]s, the `Equivalent` variant takes higher 
    /// priority, and so only `Equivalent` is used to describe the 
    /// configuration.
    Intersecting,

    /// Represents two objects not sharing any [`Point`]s while
    /// also having different directions.
    /// 
    /// Contrast this with the `Parallel` variant, where the objects 
    /// also don't share any Points but where their directions 
    /// are the same. **Note:** however, where the relative
    /// direction of the two objects is not relevant, either
    /// because they have no directions (like two Points) or
    /// because they can only have one relative direction when
    /// not sharing any Points (like a [`Line`] and a [`Plane`]),
    /// the `Parallel` variant takes priority over `Skew`, and
    /// is the only one used to describe such configurations.
    Skew
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

impl AsRef<Point> for Point {
    fn as_ref(&self) -> &Point {
        &self
    }
}

impl From<Vector> for Point {
    fn from(vector: Vector) -> Self {
        Point{x: vector.x, y: vector.y, z: vector.z}
    }
}

impl From<[f64; 3]> for Point {
    fn from(value: [f64; 3]) -> Self {
        Point{x: value[0], y: value[1], z: value[2]}
    }
}

impl From<(f64, f64, f64)> for Point {
    fn from(value: (f64, f64, f64)) -> Self {
        Point{x: value.0, y: value.1, z: value.2}
    }
}

impl TryFrom<Vec<f64>> for Point {
    type Error = String;

    fn try_from(value: Vec<f64>) -> Result<Self, Self::Error> {
        let length = value.len();

        if length == 3 {
            return Ok(Point{x: value[0], y: value[1], z: value[2]})
        } else {
            return Err(format!("A Vec can be converted into a Point only if its length is 3, but the length of the provided Vec is {}", length))
        }
    }
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

impl<'a> std::ops::Add<Point> for &'a Point {
    type Output = Point;

    fn add(self, rhs: Point) -> Self::Output {
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
        self.x += rhs.x; self.y += rhs.y; self.z += rhs.z;
    }
}

impl<'a> std::ops::AddAssign<&'a Point> for Point {
    fn add_assign(&mut self, rhs: &Point) {
        self.x += rhs.x; self.y += rhs.y; self.z += rhs.z;
    }
}

impl std::ops::AddAssign<Vector> for Point {
    fn add_assign(&mut self, rhs: Vector) {
        self.x += rhs.x; self.y += rhs.y; self.z += rhs.z;
    }
}

impl<'a> std::ops::AddAssign<&'a Vector> for Point {
    fn add_assign(&mut self, rhs: &Vector) {
        self.x += rhs.x; self.y += rhs.y; self.z += rhs.z;
    }
}

impl std::ops::Div<f64> for Point {
    type Output = Point;

    fn div(self, rhs: f64) -> Self::Output {
        Point{x: self.x / rhs, y: self.y / rhs, z: self.z / rhs}
    }
}

impl std::ops::DivAssign<f64> for Point {
    fn div_assign(&mut self, rhs: f64) {
        self.x /= rhs; self.y /= rhs; self.z /= rhs;
    }
}

impl std::ops::Index<usize> for Point {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        if index == 0 {
            return &self.x
        } else if index == 1 {
            return &self.y
        } else if index == 2 {
            return &self.z
        } else {
            panic!("Index out of bounds: Point can only be indexed with values of 0, 1, or 2, but {index} was provided.")
        }
    }
}

impl std::ops::Mul<f64> for Point {
    type Output = Point;

    fn mul(self, rhs: f64) -> Self::Output {
        Point{x: self.x * rhs, y: self.y * rhs, z: self.z * rhs}
    }
}

impl std::ops::MulAssign<f64> for Point {
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs; self.y *= rhs; self.z *= rhs;
    }
}

impl std::ops::Neg for Point {
    type Output = Point;

    fn neg(self) -> Self::Output {
        Point{ x: -self.x, y: -self.y, z: -self.z }
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

impl std::ops::Sub<Vector> for Point {
    type Output = Point;

    fn sub(self, rhs: Vector) -> Self::Output {
        Point{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl<'a> std::ops::Sub<&'a Vector> for Point {
    type Output = Point;

    fn sub(self, rhs: &Vector) -> Self::Output {
        Point{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl<'a> std::ops::Sub<Vector> for &'a Point {
    type Output = Point;

    fn sub(self, rhs: Vector) -> Self::Output {
        Point{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl std::ops::SubAssign<Point> for Point {
    fn sub_assign(&mut self, rhs: Point) {
        self.x -= rhs.x; self.y -= rhs.y; self.z -= rhs.z;
    }
}

impl<'a> std::ops::SubAssign<&'a Point> for Point {
    fn sub_assign(&mut self, rhs: &Point) {
        self.x -= rhs.x; self.y -= rhs.y; self.z -= rhs.z;
    }
}

impl std::ops::SubAssign<Vector> for Point {
    fn sub_assign(&mut self, rhs: Vector) {
        self.x -= rhs.x; self.y -= rhs.y; self.z -= rhs.z;
    }
}

impl<'a> std::ops::SubAssign<&'a Vector> for Point {
    fn sub_assign(&mut self, rhs: &Vector) {
        self.x -= rhs.x; self.y -= rhs.y; self.z -= rhs.z;
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

/// Determines whether 2 Points are exactly equal. This is the case
/// when each of their coordinates (x, y, z) are within 1e-10 of each other.
impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        return (self.x - other.x).abs() <= 1e-10 && (self.y - other.y).abs() <= 1e-10 && (self.z - other.z).abs() <= 1e-10
    }
}

/// Determines whether 2 Points are equal within a tolerance.
/// 
/// 2 Points are approximately equal when each of their coordinates
/// (x, y, z) are within the tolerance.
/// 
/// # Example
/// 
/// ```
/// use point_group::geometry::*;
/// 
/// let point = Point::new(0.0, 0.0, 0.0);
/// 
/// assert!(point.approx_eq(&Point::new(0.0, 0.0, 0.0), 1e-10));
/// 
/// assert!(!point.approx_eq(&Point::new(0.1, 0.0, -0.5), 1e-10));
/// assert!(point.approx_eq(&Point::new(0.1, 0.0, -0.5), 1.0));
/// ```
impl ApproxEq for Point {
    fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        return (self.x - other.x).abs() <= tolerance && (self.y - other.y).abs() <= tolerance && 
            (self.z - other.z).abs() <= tolerance
    }
}

impl RelativePositionTrait<Point> for Point {
    /// Computes the displacement of this Point from another
    /// Point.
    /// 
    /// **Note:** The direction of the displacement [`Vector`]
    /// is from the Point on which this method is called to
    /// the Point passed in to the method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let origin = Point::new(0.0, 0.0, 0.0);
    /// let point = Point::new(1.0, 0.0, 0.0);
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), origin.displacement(Point::new(0.0, 0.0, 0.0)));
    /// 
    /// assert_eq!(Vector::new(1.0, 0.0, 0.0), origin.displacement(point));
    /// assert_eq!(Vector::new(-1.0, 0.0, 0.0), point.displacement(origin));
    /// 
    /// assert_eq!(point.displacement(origin), - origin.displacement(point));
    /// ```
    fn displacement(self, other: Point) -> Vector {
        Vector::from_two_points(self, other)
    }

    /// Computes the distance between 2 Points.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let origin = Point::new(0.0, 0.0, 0.0);
    /// 
    /// assert_eq!(0.0, origin.distance(Point::new(0.0, 0.0, 0.0)));
    /// assert_eq!(1.0, origin.distance(Point::new(1.0, 0.0, 0.0)));
    /// assert_eq!(5.0, origin.distance(Point::new(0.0, -3.0, 4.0)));
    /// 
    /// // The distance is equal regardless of which Point the method is calledo on.
    /// let point = Point::new(-1.5, 10.0, 0.1);
    /// assert_eq!(origin.distance(point), point.distance(origin));
    /// ```
    fn distance(self, other: Point) -> f64 {
        self.displacement(other).magnitude()
    }

    /// Determines the relative position of a Point and a Point. 
    /// 
    /// Two Points can adopt the following [`RelativePosition`]s in 3D space:
    ///  1. Equivalent
    ///  2. Parallel
    /// 
    /// Two Points are Equivalent when they are exactly equal.
    /// 
    /// Two Points are Parallel when they are not exactly equal.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let origin = Point::new(0.0, 0.0, 0.0);
    /// 
    /// assert_eq!(RelativePosition::Equivalent, origin.relative_position(Point::new(0.0, 0.0, 0.0)));
    /// assert_eq!(RelativePosition::Parallel, origin.relative_position(Point::new(0.0, -3.0, 4.0)));
    /// ```
    fn relative_position(self, other: Point) -> RelativePosition {
        if self == other {
            return RelativePosition::Equivalent
        } else {
            return RelativePosition::Parallel
        }
    }
}

impl RelativePositionTrait<Line> for Point {
    /// Computes the displacement of this Point from a [`Line`].
    /// 
    /// **Note:** The displacement [`Vector`] is in the direction
    /// from the Point to the Line.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), Point::new(0.5, 0.0, 0.0).displacement(line));
    /// assert_eq!(Vector::new(0.0, -1.0, 0.0), Point::new(0.0, 1.0, 0.0).displacement(line));
    /// assert_eq!(Vector::new(0.0, 1.0, 0.0), Point::new(0.0, -1.0, 0.0).displacement(line));
    /// ```
    fn displacement(self, line: Line) -> Vector {
        let cross = line.vector.cross_product(Vector::from_two_points(self, line.point));

        let distance = cross.magnitude() / line.vector.magnitude();

        if distance < 1e-10 {
            return Vector::new(0.0, 0.0, 0.0)
        }

        let vector = cross.cross_product(line.vector);

        return vector * distance / vector.magnitude()
    }

    /// Computes the distance of this Point from a [`Line`]. The following
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
    /// assert_eq!(0.0, Point::new(0.5, 0.0, 0.0).distance(line));
    /// assert_eq!(1.0, Point::new(0.0, 1.0, 0.0).distance(line));
    /// ```
    fn distance(self, line: Line) -> f64 {
        line.vector.cross_product(Vector::from_two_points(self, line.point)).magnitude() / line.vector.magnitude()
    }

    /// Determines the relative position of a Point and a [`Line`].
    /// 
    /// A Point and a Line can adopt the following [`RelativePosition`]s
    /// in 3D space:
    ///  1. Equivalent
    ///  2. Parallel
    /// 
    /// A Point and a Line are Equivalent when the Point lies on the
    /// Line.
    /// 
    /// They are Parallel when the Point does not lie on the Line.
    /// 
    /// # Example
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert_eq!(RelativePosition::Equivalent, Point::new(0.5, 0.0, 0.0).relative_position(line));
    /// assert_eq!(RelativePosition::Parallel, Point::new(0.0, 1.0, 0.0).relative_position(line));
    /// ```
    fn relative_position(self, line: Line) -> RelativePosition {
        if self.is_on_line(line) {
            return RelativePosition::Equivalent
        } else {
            return RelativePosition::Parallel
        }
    }
}

impl RelativePositionTrait<Plane> for Point {
    /// Computes the displacement of this Point from a [`Plane`]. This 
    /// is done by multiplying the normal [`Vector`] by the negative
    /// distance of the Point from the Plane divided by the magnitude
    /// of the normal.
    /// 
    /// **Note:** The direction of the returned displacement Vector
    /// is from the Point towards the Plane. fhis is the opposite 
    /// displacement to the one obtained from the `Plane::displacement`
    /// method, i.e.: 
    /// `point.displacement(plane) == - plane.displacement(point)`
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), Point::new(0.0, 5.0, -1.0).displacement(plane));
    /// 
    /// // When the Vector from the Point to the Plane has the opposite direction to the
    /// // normal Vector, the displacement Vector will be a negative k-multiple.
    /// assert_eq!(Vector::new(-1.0, 0.0, 0.0), Point::new(1.0, 5.0, -1.0).displacement(plane));
    /// 
    /// // When the Vector from the Point to the Plane has the same direction to the
    /// // normal Vector, the displacement Vector will be a positive k-multiple.
    /// assert_eq!(Vector::new(1.0, 0.0, 0.0), Point::new(-1.0, 5.0, -1.0).displacement(plane));
    /// 
    /// // Point::displacement() and Plane::displacement() are opposites.
    /// let point = Point::new(3.0, 0.0, 0.0);
    /// assert_eq!(point.displacement(plane), - plane.displacement(point));
    /// ```
    fn displacement(self, plane: Plane) -> Vector {
        let magnitude = plane.normal.magnitude();
        let distance = Vector::from_two_points(plane.point, self).dot_product(plane.normal) / magnitude;

        return plane.normal * distance * -1.0 / magnitude
    }

    /// Computes the distance of this Point from a [`Plane`]. The following
    /// equation is used:
    ///     [(x - y) . n] / |n|
    /// where x is this Point, y is the point defining the Plane, and n is
    /// the normal of the Plane.
    /// 
    /// **Note:** The computed distance is positive if the normal [`Vector`] points
    /// to the same side of the Plane as the one on which this Point is located
    /// and vice versa. In other words, if the k multiple between the normal 
    /// vector from the Plane to the Point and the normal vector defining the 
    /// Plane is positive (they are in the same direction), the distance will
    /// be positive. If the k multiple is negative (they are in the opposite
    /// directions), the distance will be negative.
    /// 
    /// **Note to the Note:** This positive/negative value caused by the 
    /// relative positions of the Point and the Plane is the same as when
    /// this method is called on Plane: 
    /// `point.distance(plane) == plane.distance(point)`
    /// 
    /// To obtain the absolute distance, simply call the abs() method: 
    /// `point.distance(plane).abs()`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(0.0, Point::new(0.0, 5.0, -1.0).distance(plane));
    /// 
    /// // The normal vector is parallel to the x-axis, and the Point is displaced
    /// // to the right along the x-axis (towards positive infinity), so the computed
    /// // distance is positive.
    /// assert_eq!(1.0, Point::new(1.0, 5.0, -1.0).distance(plane));
    /// 
    /// // The normal vector is parallel to the x-axis, and the Point is displaced
    /// // to the left along the x-axis (towards negative infinity), so the computed
    /// // distance is negative.
    /// assert_eq!(-1.0, Point::new(-1.0, 5.0, -1.0).distance(plane));
    /// 
    /// // The distance from a Point to a Plane and vice versa are equivalent.
    /// let point = Point::new(2.0, 0.0, 0.0);
    /// assert_eq!(point.distance(plane), plane.distance(point));
    /// ```
    fn distance(self, plane: Plane) -> f64 {
        Vector::from_two_points(plane.point, self).dot_product(plane.normal) / plane.normal.magnitude()
    }

    /// Determines the relative position of a Point and a [`Plane`].
    /// 
    /// A Point and a Plane can adopt the following [`RelativePosition`]s 
    /// in 3D space:
    ///  1. Equivalent
    ///  2. Parallel
    /// 
    /// A Point and a Plane are Equivalent when the Point lies on the
    /// Plane.
    /// 
    /// They are Parallel when the Point does not lie on the Plane.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(RelativePosition::Equivalent, Point::new(0.0, 5.0, -1.0).relative_position(plane));
    /// assert_eq!(RelativePosition::Parallel, Point::new(1.0, 5.0, -1.0).relative_position(plane));
    /// ```
    fn relative_position(self, plane: Plane) -> RelativePosition {
        if self.is_on_plane(plane) {
            return RelativePosition::Equivalent
        } else {
            return RelativePosition::Parallel
        }
    }
}

impl VectorAlgebra<Point> for Point {
    type Output = Point;

    fn cross_product(self, other: Point) -> Self::Output {
        let (x, y, z) = _cross_product(self.x, self.y, self.z, other.x, other.y, other.z);

        return Point{x, y, z}
    }

    fn dot_product(self, other: Point) -> f64 {
        _dot_product(self.x, self.y, self.z, other.x, other.y, other.z)
    }
}

impl VectorAlgebra<Vector> for Point {
    type Output = Vector;

    fn cross_product(self, other: Vector) -> Self::Output {
        let (x, y, z) = _cross_product(self.x, self.y, self.z, other.x, other.y, other.z);

        return Vector{x, y, z}
    }

    fn dot_product(self, other: Vector) -> f64 {
        _dot_product(self.x, self.y, self.z, other.x, other.y, other.z)
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
        self.distance(plane).abs() < tolerance
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
        let displacement = self.displacement(plane);

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

impl std::ops::Add<Point> for Vector {
    type Output = Vector;

    fn add(self, rhs: Point) -> Self::Output {
        Vector{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl<'a> std::ops::Add<&'a Point> for Vector {
    type Output = Vector;

    fn add(self, rhs: &Point) -> Self::Output {
        Vector{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl<'a> std::ops::Add<Point> for &'a Vector {
    type Output = Vector;

    fn add(self, rhs: Point) -> Self::Output {
        Vector{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl std::ops::Add<Vector> for Vector {
    type Output = Vector;

    fn add(self, rhs: Vector) -> Self::Output {
        Vector{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl<'a> std::ops::Add<&'a Vector> for Vector {
    type Output = Vector;

    fn add(self, rhs: &Vector) -> Self::Output {
        Vector{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl<'a> std::ops::Add<Vector> for &'a Vector {
    type Output = Vector;

    fn add(self, rhs: Vector) -> Self::Output {
        Vector{x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z}
    }
}

impl std::ops::AddAssign<Point> for Vector {
    fn add_assign(&mut self, rhs: Point) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<'a> std::ops::AddAssign<&'a Point> for Vector {
    fn add_assign(&mut self, rhs: &Point) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl std::ops::AddAssign<Vector> for Vector {
    fn add_assign(&mut self, rhs: Vector) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<'a> std::ops::AddAssign<&'a Vector> for Vector {
    fn add_assign(&mut self, rhs: &Vector) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

/// Determines whether 2 Vectors are equal within a tolerance.
/// 
/// 2 Vectors are approximately equal when the components along
/// each of the axes (x, y, z) are within a tolerance.
/// 
/// # Example
/// 
/// ```
/// use point_group::geometry::*;
/// 
/// let vector = Vector::new(0.0, 0.0, 0.0);
/// 
/// assert!(vector.approx_eq(&Vector::new(0.0, 0.0, 0.0), 1e-10));
/// 
/// assert!(!vector.approx_eq(&Vector::new(0.1, 0.0, -0.5), 1e-10));
/// assert!(vector.approx_eq(&Vector::new(0.1, 0.0, -0.5), 1.0));
/// ```
impl ApproxEq for Vector {
    fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        return (self.x - other.x).abs() <= tolerance && (self.y - other.y).abs() <= tolerance && 
            (self.z - other.z).abs() <= tolerance
    }
}

impl AsRef<Vector> for Vector {
    fn as_ref(&self) -> &Vector {
        &self
    }
}

impl std::ops::Div<f64> for Vector {
    type Output = Vector;

    fn div(self, rhs: f64) -> Self::Output {
        Vector{x: self.x / rhs, y: self.y / rhs, z: self.z / rhs}
    }
}

impl std::ops::DivAssign<f64> for Vector {
    fn div_assign(&mut self, rhs: f64) {
        self.x /= rhs; self.y /= rhs; self.z /= rhs;
    }
}

impl From<Point> for Vector {
    fn from(point: Point) -> Self {
        Vector{x: point.x, y: point.y, z: point.z}
    }
}

impl From<[f64; 3]> for Vector {
    fn from(value: [f64; 3]) -> Self {
        Vector{x: value[0], y: value[1], z: value[2]}
    }
}

impl From<(f64, f64, f64)> for Vector {
    fn from(value: (f64, f64, f64)) -> Self {
        Vector{x: value.0, y: value.1, z: value.2}
    }
}

impl TryFrom<Vec<f64>> for Vector {
    type Error = String;

    fn try_from(value: Vec<f64>) -> Result<Self, Self::Error> {
        let length = value.len();

        if length == 3 {
            return Ok(Vector{x: value[0], y: value[1], z: value[2]})
        } else {
            return Err(format!("A Vec can be converted into a Vector only if its length is 3, but the length of the provided Vec is {}", length))
        }
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
        self.x *= rhs; self.y *= rhs; self.z *= rhs;
    }
}

impl std::ops::Neg for Vector {
    type Output = Vector;

    fn neg(self) -> Self::Output {
        Vector { x: -self.x, y: -self.y, z: -self.z }
    }
}

/// Determines whether 2 Vectors are exactly equal. This is the case
/// when both their components along each of the axes (x, y, z) are 
/// within 1e-10 of each other.
impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        return (self.x - other.x).abs() <= 1e-10 && (self.y - other.y).abs() <= 1e-10 && (self.z - other.z).abs() <= 1e-10
    }
}

/// Determines the relative size of 2 Vectors by their magnitudes.
/// 
/// # Example
/// 
/// ```
/// use point_group::geometry::*;
/// 
/// let vector1 = Vector::new(1.0, 0.0, 0.0);
/// let vector2 = Vector::new(-1.0, 0.0, 0.0);
/// let vector3 = Vector::new(0.0, 0.0, 10.0);
/// 
/// assert!(vector3 > vector1);
/// assert!(vector2 < vector3);
/// 
/// assert!(vector1 <= vector2);
/// assert!(vector1 >= vector2);
/// ```
impl PartialOrd for Vector {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let self_mag = self.magnitude();
        let other_mag = other.magnitude();

        if (self_mag - other_mag).abs() < 1e-10 {
            return Some(std::cmp::Ordering::Equal)
        } else if self_mag > other_mag {
            return Some(std::cmp::Ordering::Greater)
        } else if self_mag < other_mag {
            return Some(std::cmp::Ordering::Less)
        } else {
            None
        }
    }
}

impl std::ops::Index<usize> for Vector {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        if index == 0 {
            return &self.x
        } else if index == 1 {
            return &self.y
        } else if index == 2 {
            return &self.z
        } else {
            panic!("Index out of bounds: Point can only be indexed with values of 0, 1, or 2, but {index} was provided.")
        }
    }
}

impl std::ops::Sub<Point> for Vector {
    type Output = Vector;

    fn sub(self, rhs: Point) -> Self::Output {
        Vector{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl<'a> std::ops::Sub<&'a Point> for Vector {
    type Output = Vector;

    fn sub(self, rhs: &Point) -> Self::Output {
        Vector{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl<'a> std::ops::Sub<Point> for &'a Vector {
    type Output = Vector;

    fn sub(self, rhs: Point) -> Self::Output {
        Vector{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl std::ops::Sub<Vector> for Vector {
    type Output = Vector;

    fn sub(self, rhs: Vector) -> Self::Output {
        Vector{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl<'a> std::ops::Sub<&'a Vector> for Vector {
    type Output = Vector;

    fn sub(self, rhs: &Vector) -> Self::Output {
        Vector{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl<'a> std::ops::Sub<Vector> for &'a Vector {
    type Output = Vector;

    fn sub(self, rhs: Vector) -> Self::Output {
        Vector{x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z}
    }
}

impl std::ops::SubAssign<Point> for Vector {
    fn sub_assign(&mut self, rhs: Point) {
        self.x -= rhs.x; self.y -= rhs.y; self.z -= rhs.z;
    }
}

impl<'a> std::ops::SubAssign<&'a Point> for Vector {
    fn sub_assign(&mut self, rhs: &Point) {
        self.x -= rhs.x; self.y -= rhs.y; self.z -= rhs.z;
    }
}

impl std::ops::SubAssign<Vector> for Vector {
    fn sub_assign(&mut self, rhs: Vector) {
        self.x -= rhs.x; self.y -= rhs.y; self.z -= rhs.z;
    }
}

impl<'a> std::ops::SubAssign<&'a Vector> for Vector {
    fn sub_assign(&mut self, rhs: &Vector) {
        self.x -= rhs.x; self.y -= rhs.y; self.z -= rhs.z;
    }
}

impl std::iter::Sum<Vector> for Vector {
    fn sum<I: Iterator<Item = Vector>>(iter: I) -> Self {
        iter.fold(Vector{x: 0.0, y: 0.0, z: 0.0}, |acc, x| acc + x)
    }
}

impl<'a> std::iter::Sum<&'a Vector> for Vector {
    fn sum<I: Iterator<Item = &'a Vector>>(iter: I) -> Self {
        iter.fold(Vector{x: 0.0, y: 0.0, z: 0.0}, |acc, x| acc + *x)
    }
}

impl VectorAlgebra<Point> for Vector {
    type Output = Vector;

    fn cross_product(self, other: Point) -> Self::Output {
        let (x, y, z) = _cross_product(self.x, self.y, self.z, other.x, other.y, other.z);

        return Vector{x, y, z}
    }

    fn dot_product(self, other: Point) -> f64 {
        _dot_product(self.x, self.y, self.z, other.x, other.y, other.z)
    }
}

impl VectorAlgebra<Vector> for Vector {
    type Output = Vector;

    fn cross_product(self, other: Vector) -> Self::Output {
        let (x, y, z) = _cross_product(self.x, self.y, self.z, other.x, other.y, other.z);

        return Vector{x, y, z}
    }

    fn dot_product(self, other: Vector) -> f64 {
        _dot_product(self.x, self.y, self.z, other.x, other.y, other.z)
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

    /// Returns a normalised (i.e. divided by its length) copy of this Vector.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let normalised_vector = Vector::new(1.0, 0.0, 0.0);
    /// let k_multiple = Vector::new(10.0, 0.0, 0.0);
    /// let random_vector = Vector::new(0.0, 4.0, -3.0);
    /// 
    /// assert_eq!(Vector::new(1.0, 0.0, 0.0), normalised_vector.normalise());
    /// assert_eq!(Vector::new(1.0, 0.0, 0.0), k_multiple.normalise());
    /// assert_eq!(Vector::new(0.0, 0.8, -0.6), random_vector.normalise())
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
/// 
/// In this module, this struct is used solely for 3D rotation 
/// computations, and so its mathematical implementation is 
/// incomplete.
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

impl AsRef<Quaternion> for Quaternion {
    fn as_ref(&self) -> &Quaternion {
        &self
    }
}

impl AsRef<Vector> for Quaternion {
    fn as_ref(&self) -> &Vector {
        &self.vector
    }
}

impl From<[f64; 4]> for Quaternion {
    fn from(value: [f64; 4]) -> Self {
        Quaternion { q: value[0], vector: Vector { x: value[1], y: value[2], z: value[3] } }
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

impl AsRef<Line> for Line {
    fn as_ref(&self) -> &Line {
        &self
    }
}

impl AsRef<Point> for Line {
    fn as_ref(&self) -> &Point {
        &self.point
    }
}

impl AsRef<Vector> for Line {
    fn as_ref(&self) -> &Vector {
        &self.vector
    }
}

impl PartialEq for Line {
    fn eq(&self, other: &Self) -> bool {
        return self.vector.is_k_multiple(other.vector) && other.point.is_on_line(*self)
    }
}

/// Used for determining the relative positions of a Line and a [`Point`].
impl RelativePositionTrait<Point> for Line {
    /// Computes the displacement of this Line from a [`Point`].
    /// 
    /// **Note:** The displacement [`Vector`] is in the direction
    /// from the Line to the Point.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), line.displacement(Point::new(0.5, 0.0, 0.0)));
    /// assert_eq!(Vector::new(0.0, 1.0, 0.0), line.displacement(Point::new(0.0, 1.0, 0.0)));
    /// assert_eq!(Vector::new(0.0, -1.0, 0.0), line.displacement(Point::new(0.0, -1.0, 0.0)));
    /// ```
    fn displacement(self, point: Point) -> Vector {
        let cross = self.vector.cross_product(Vector::from_two_points(self.point, point));

        let distance = cross.magnitude() / self.vector.magnitude();

        if distance < 1e-10 {
            return Vector::new(0.0, 0.0, 0.0)
        }

        let vector = cross.cross_product(self.vector);

        return vector * distance / vector.magnitude()
    }

    /// Computes the distance of this Line from a [`Point`]. 
    /// 
    /// This method is just a wrapper around the `Point::distance()`
    /// method.
    /// 
    /// # Example
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert_eq!(0.0, line.distance(Point::new(0.5, 0.0, 0.0)));
    /// assert_eq!(1.0, line.distance(Point::new(0.0, 1.0, 0.0)));
    /// ```
    fn distance(self, point: Point) -> f64 {
        point.distance(self)
    }

    /// Determines the relative position of a Line and a [`Point`].
    /// 
    /// A Line and a Point can adopt the following [`RelativePosition`]s
    /// in 3D space:
    ///  1. Equivalent
    ///  2. Parallel
    /// 
    /// A Line and a Point are Equivalent when the Point lies on the
    /// Line.
    /// 
    /// They are Parallel when the Point does not lie on the Line.
    /// 
    /// # Example
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let line = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// assert_eq!(RelativePosition::Equivalent, Point::new(0.5, 0.0, 0.0).relative_position(line));
    /// assert_eq!(RelativePosition::Parallel, Point::new(0.0, 1.0, 0.0).relative_position(line));
    /// ```
    fn relative_position(self, point: Point) -> RelativePosition {
        point.relative_position(self)
    }
}

/// Used for determining the relative positions of 2 Lines
impl RelativePositionTrait<Line> for Line {
    /// Computes the displacement between 2 Lines. A 0 [`Vector`] is returned
    /// when the Lines' [`RelativePosition`] is Equivalent or Intersecting.
    /// 
    /// **Note:** The direction of the displacement Vector is from the Line
    /// this method is called on to the Line passed in to this method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// let y = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0));
    /// let parallel_x = Line::new(Point::new(0.0, 0.1, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// let parallel_y = Line::new(Point::new(0.0, 0.0, -0.1), Vector::new(0.0, 1.0, 0.0));
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), x.displacement(x), "xx");
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), x.displacement(y), "yy");
    /// 
    /// assert_eq!(Vector::new(0.0, 0.1, 0.0), x.displacement(parallel_x));
    /// assert_eq!(Vector::new(0.0, -0.1, 0.0), parallel_x.displacement(x));
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, -0.1), x.displacement(parallel_y));
    /// ```
    fn displacement(self, other: Line) -> Vector {
        match self.relative_position(other) {
            RelativePosition::Equivalent | RelativePosition::Intersecting => return Vector::new(0.0, 0.0, 0.0),
            RelativePosition::Parallel => return self.displacement(other.point),
            RelativePosition::Skew => {
                let cross = self.vector.cross_product(other.vector);
                let magnitude = cross.magnitude();
                let distance = Vector::from_two_points(self.point, other.point).dot_product(cross) / magnitude;
                return cross * distance / magnitude
            }
        }
    }

    /// Computes the distance between 2 Lines. 0.0 is returned when
    /// the Lines' [`RelativePosition`] is Equivalent or Intersecting,
    /// otherwise the distance is computed.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// let y = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0));
    /// let parallel_x = Line::new(Point::new(0.0, 0.1, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// let parallel_y = Line::new(Point::new(0.0, 0.0, -0.5), Vector::new(0.0, 1.0, 0.0));
    /// 
    /// assert_eq!(0.0, x.distance(x));
    /// assert_eq!(0.0, x.distance(y));
    /// assert_eq!(0.1, x.distance(parallel_x));
    /// assert_eq!(0.5, x.distance(parallel_y));
    /// ```
    fn distance(self, other: Line) -> f64 {
        match self.relative_position(other) {
            RelativePosition::Equivalent | RelativePosition::Intersecting => 0.0,
            RelativePosition::Parallel => self.distance(other.point),
            RelativePosition::Skew => {
                let cross = self.vector.cross_product(other.vector);
                return Vector::from_two_points(self.point, other.point).dot_product(cross).abs() / cross.magnitude()
            }
        }
    }

    /// Determines the relative position of two Lines in space.
    /// 
    /// Two Lines can adopt these [`RelativePosition`]s in 3D space:
    ///  1. Equivalent
    ///  2. Parallel
    ///  3. Intersecting
    ///  4. Skew
    /// 
    /// Twp Lines are Equivalent when their directional [`Vector`]s are
    /// k multiples of one another and they share the [`Point`]s that 
    /// define them.
    /// 
    /// They are Parallel when their directional [`Vector`]s are
    /// k multiples of one another but they do not share the [`Point`]s  
    /// that define them.
    /// 
    /// They are Intersecting  when their directional [`Vector`]s are
    /// not k multiples of one another and their parametric equations
    /// are equal.
    /// 
    /// They are Skew  when their directional [`Vector`]s are not
    /// k multiples of one another and their parametric equations
    /// are not equal.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// assert_eq!(RelativePosition::Equivalent, x.relative_position(x));
    /// 
    /// let parallel_x = Line::new(Point::new(0.0, 0.1, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// assert_eq!(RelativePosition::Parallel, x.relative_position(parallel_x));
    /// 
    /// let y = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0));
    /// assert_eq!(RelativePosition::Intersecting, x.relative_position(y));
    /// 
    /// let parallel_y = Line::new(Point::new(0.0, 0.0, -0.1), Vector::new(0.0, 1.0, 0.0));
    /// assert_eq!(RelativePosition::Skew, x.relative_position(parallel_y));
    /// ```
    fn relative_position(self, other: Line) -> RelativePosition {
        if self.vector.is_k_multiple(other.vector) {
            if self.has_point(other.point) {
                return RelativePosition::Equivalent
            } else {
                return RelativePosition::Parallel
            }
        } else {
            let self_vec = [self.vector.x.abs(), self.vector.y.abs(), self.vector.z.abs()];
            let other_vec = [other.vector.x.abs(), other.vector.y.abs(), other.vector.z.abs()];

            for (i, (self_val, other_val)) in self_vec.iter().zip(other_vec.iter()).enumerate() {
                if self_val < &1e-10 && other_val < &1e-10 {
                    if (self.point[i] - other.point[i]).abs() < 1e-10 {
                        return RelativePosition::Intersecting
                    } else {
                        return RelativePosition::Skew
                    }
                }
            }

            for (i, (self_val, other_val)) in self_vec[..2].iter().zip(other_vec[..2].iter()).enumerate() {
                let j = i + 1;
                let k = if i == 0 { i + 2} else { i - 1 };

                let (mu, lamda, pos) = if *self_val < 1e-10 && *other_val > 1e-10 {
                    let mu = (self.point[i] - other.point[i]) / other.vector[i];
                    let (lamda, pos) = if self_vec[j] > 1e-10 {
                        ((other.point[j] - self.point[j] + mu * other.vector[j]) / self.vector[j], k)
                    } else {
                        ((other.point[k] - self.point[k] + mu * other.vector[k]) / self.vector[k], j)
                    };
                    (mu, lamda, pos)
                } else if *self_val > 1e-10 && *other_val < 1e-10 {
                    let lamda = (other.point[i] - self.point[i]) / self.vector[i];
                    let (mu, pos) = if other_vec[j] > 1e-10 {
                        ((- other.point[j] + self.point[j] + lamda * self.vector[j]) / other.vector[j], k)
                    } else {
                        ((other.point[k] - self.point[k] + lamda * self.vector[k]) / other.vector[k], j)
                    };
                    (mu, lamda, pos)
                } else {
                    (std::f64::NAN, std::f64::NAN, 0)
                };

                if !mu.is_nan() && !lamda.is_nan() {
                    if ((self.point[pos] + self.vector[pos] * lamda) - (other.point[pos] + other.vector[pos] * mu)).abs() < 1e-10 {
                        return RelativePosition::Intersecting
                    } else {
                        return RelativePosition::Skew
                    }
                }
            }

            let mu = (self.vector.y * (self.point.x - other.point.x) + self.vector.x * (other.point.y - self.point.y)) / 
                (self.vector.y * other.vector.x - self.vector.x * other.vector.y);
            let lamda = (other.point.y + other.vector.y * mu - self.point.y) / self.vector.y;

            if self.point.z + self.vector.z * lamda == other.point.z + other.vector.z * mu {
                return RelativePosition::Intersecting
            } else {
                return RelativePosition::Skew
            }
        }
    }
}

/// Used for determining the relative positions of a Line and a [`Plane`].
impl RelativePositionTrait<Plane> for Line {
    /// Computes the displacement between a Line and a [`Plane`]. A zero [`Vector`]
    /// is returned if their [`RelativePosition`] is Equivalent or Parallel.
    /// 
    /// **Note:** The computed displacement Vector is in the direction from the
    /// Line towards the Plane.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// // The x axis belongs to the xy Plane
    /// let xy = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), x.displacement(xy));
    /// 
    /// // The x axis is parallel to the Plane parallel to xy
    /// let parallel_xy = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(Vector::new(0.0, 0.0, 2.0), x.displacement(parallel_xy));
    /// 
    /// // The x axis intersects with the yz Plane in exactly one Point
    /// let yz = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), x.displacement(yz));
    /// ```
    fn displacement(self, plane: Plane) -> Vector {
        match self.relative_position(plane) {
            RelativePosition::Parallel => return self.point.displacement(plane),
            _ => return Vector::new(0.0, 0.0, 0.0)
        }
    }

    /// Computes the distance between a Line and a [`Plane`]. 0.0 is returned
    /// when their [`RelativePosition`] is Equivalent of Parallel.
    /// 
    /// **Note:** Unlike the `Point::distance(plane: Plane)` method, the 
    /// distance computed by this method is always positive.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// // The x axis belongs to the xy Plane
    /// let xy = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(0.0, x.distance(xy));
    /// 
    /// // The x axis is parallel to the Plane parallel to xy
    /// let parallel_xy = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(2.0, x.distance(parallel_xy));
    /// 
    /// // The x axis intersects with the yz Plane in exactly one Point
    /// let yz = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(0.0, x.distance(yz));
    /// ```
    fn distance(self, plane: Plane) -> f64 {
        match self.relative_position(plane) {
            RelativePosition::Parallel => return self.point.distance(plane).abs(),
            _ => return 0.0
        }
    }

    /// Determines the relative position of a Line and a [`Plane`].
    /// 
    /// A Line and a Plane can adopt these [`RelativePosition`]s in 3D space:
    ///  1. Equivalent
    ///  2. Parallel
    ///  3. Intersecting
    /// 
    /// A Line and a Plane are Equivalent if the Line belongs to the Plane,
    /// that is the Line's directional [`Vector`] and the Plane's normal Vector
    /// are perpendicular and the [`Point`] defining the Line lies on the Plane.
    /// 
    /// A Line and a Plane are Parallel if they do not share any Points, i.e.
    /// if their Vectors are perpendicular but the Line's Point does not lie on
    /// the Plane.
    /// 
    /// A Line and a Plane are Intersecting if they share exactly one Point. 
    /// That is the case when their Vectors are not perpendicular.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// // The x axis belongs to the xy Plane
    /// let xy = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(RelativePosition::Equivalent, x.relative_position(xy));
    /// 
    /// // The x axis is parallel to the Plane parallel to xy
    /// let parallel_xy = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(RelativePosition::Parallel, x.relative_position(parallel_xy));
    /// 
    /// // The x axis intersects with the yz Plane in exactly one Point
    /// let yz = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(RelativePosition::Intersecting, x.relative_position(yz));
    /// ```
    fn relative_position(self, plane: Plane) -> RelativePosition {
        if self.vector.dot_product(plane.normal).abs() < 1e-10 {
            if self.has_point(plane.point) {
                return RelativePosition::Equivalent
            } else {
                return RelativePosition::Parallel
            }
        } else {
            return RelativePosition::Intersecting
        }
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

    /// Determines wheter two Lines are parallel. Two Lines are
    /// parallel when their direction [`Vector`]s are k multiples
    /// of each other.
    /// 
    /// Further, if it is known that the two Lines share a common
    /// [`Point`], this method can also be used to check their 
    /// equivalence.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// // Parallel Lines return true as long as they are within tolerance
    /// let parallel = Line::new(Point::new(0.0, -3.0, 7.5), Vector::new(1.0, 1e-5, 0.0));
    /// assert!(x.fast_approx_eq(parallel, 1e-3));
    /// assert!(!x.approx_eq(&parallel, 1e-3));
    /// 
    /// // Therefore, if it is known they share a Point, it is certain they are equivalent
    /// let almost_x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 1e-5, 0.0));
    /// assert!(x.fast_approx_eq(almost_x, 1e-3));
    /// assert!(x.approx_eq(&almost_x, 1e-3));
    /// 
    /// // Intercepting Lines return false
    /// let y = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0));
    /// assert!(!x.fast_approx_eq(y, 1e-3));
    /// ```
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

    /// Creates a copy of this Line but with its direction [`Vector`] normalised.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let already_normalised = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// let non_normalised = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(5.0, 0.0, 0.0));
    /// 
    /// assert_eq!(already_normalised, already_normalised.normalise());
    /// assert_eq!(already_normalised, non_normalised.normalise());
    /// ```
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

impl AsRef<Plane> for Plane {
    fn as_ref(&self) -> &Plane {
        &self
    }
}

impl AsRef<Point> for Plane {
    fn as_ref(&self) -> &Point {
        &self.point
    }
}

impl AsRef<Vector> for Plane {
    fn as_ref(&self) -> &Vector {
        &self.normal
    }
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

/// Used for determining the realtive position of a Plane and a [`Point`].
impl RelativePositionTrait<Point> for Plane {
    /// Computes the displacement of this Plane from a [`Point`]. This 
    /// is done by multiplying the normal [`Vector`] by the distance
    /// of the Point from the Plane divided by the magnitude of the 
    /// normal.
    /// 
    /// **Note:** The direction of the returned displacement Vector
    /// is from the Plane towards the Point. This is the opposite 
    /// displacement to the one obtained from the `Point::displacement`
    /// method, i.e.: 
    /// `plane.displacement(point) == - point.displacement(plane)`
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), plane.displacement(Point::new(0.0, 5.0, -1.0)));
    /// 
    /// // When the Vector from the Plane to the Point has the same direction as the
    /// // normal Vector, the displacement Vector will be a positive k-multiple.
    /// assert_eq!(Vector::new(1.0, 0.0, 0.0), plane.displacement(Point::new(1.0, 5.0, -1.0)));
    /// 
    /// // When the Vector from the Plane to the Point has the opposite direction to the
    /// // normal Vector, the displacement Vector will be a negative k-multiple.
    /// assert_eq!(Vector::new(-2.0, 0.0, 0.0), plane.displacement(Point::new(-2.0, 5.0, -1.0)));
    /// 
    /// // Plane::displacement() and Point::displacement() are opposites.
    /// let point = Point::new(3.0, 0.0, 0.0);
    /// assert_eq!(plane.displacement(point), - point.displacement(plane));
    /// ```
    fn displacement(self, point: Point) -> Vector {
        let distance = self.distance(point);
        return self.normal * distance / self.normal.magnitude()
    }

    /// Computes the distance of this Plane from a [`Point`]. The following
    /// equation is used:
    ///     [(x - y) . n] / |n|
    /// where x is the Point, y is the Point defining this Plane, and n is
    /// the normal of the Plane.
    /// 
    /// **Note:** The computed distance is positive if the normal [`Vector`] points
    /// to the same side of the Plane as the one on which this Point is located
    /// and vice versa. In other words, if the k multiple between the normal 
    /// vector from the Plane to the Point and the normal vector defining the 
    /// Plane is positive (they are in the same direction), the distance will
    /// be positive. If the k multiple is negative (they are in the opposite
    /// directions), the distance will be negative.
    /// 
    /// **Note to the Note:** This positive/negative value caused by the 
    /// relative positions of the Point and the Plane is the same as when
    /// this method is called on Point: 
    /// `plane.distance(point) == point.distance(plane)`
    /// 
    /// To obtain the absolute distance, simply call the abs() method: 
    /// `plane.distance(point).abs()`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(0.0, plane.distance(Point::new(0.0, 5.0, -1.0)));
    /// 
    /// // The normal vector is parallel to the x-axis, and the Point is displaced
    /// // to the right along the x-axis (towards positive infinity), so the computed
    /// // distance is positive.
    /// assert_eq!(1.0, plane.distance(Point::new(1.0, 5.0, -1.0)));
    /// 
    /// // The normal vector is parallel to the x-axis, and the Point is displaced
    /// // to the left along the x-axis (towards negative infinity), so the computed
    /// // distance is negative.
    /// assert_eq!(-1.0, plane.distance(Point::new(-1.0, 5.0, -1.0)));
    /// 
    /// // The distance from a Plane to a Point and vice versa are equivalent.
    /// let point = Point::new(2.0, 0.0, 0.0);
    /// assert_eq!(plane.distance(point), point.distance(plane));
    /// ```
    fn distance(self, point: Point) -> f64 {
        Vector::from_two_points(self.point, point).dot_product(self.normal) / self.normal.magnitude()
    }

    /// Determines the relative position of a Plane and a [`Point`].
    /// 
    /// A Plane and a Point can adopt the following [`RelativePosition`]s 
    /// in 3D space:
    ///  1. Equivalent
    ///  2. Parallel
    /// 
    /// A Plane and a Point are Equivalent when the Point lies on the
    /// Plane.
    /// 
    /// They are Parallel when the Point does not lie on the Plane.
    /// 
    /// Note: this method is a wrapper around `Point::relative_position(plane: Plane)`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let plane = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// 
    /// assert_eq!(RelativePosition::Equivalent, plane.relative_position(Point::new(0.0, 5.0, -1.0)));
    /// assert_eq!(RelativePosition::Parallel, plane.relative_position(Point::new(1.0, 5.0, -1.0)));
    /// ```
    fn relative_position(self, point: Point) -> RelativePosition {
        point.relative_position(self)
    }
}

/// Used for determining the relative position of a Plane and a [`Line`].
impl RelativePositionTrait<Line> for Plane {
    /// Computes the displacement between a Plane and a [`Line`]. A zero [`Vector`]
    /// is returned if their [`RelativePosition`] is Equivalent or Parallel.
    /// 
    /// **Note:** The computed displacement Vector is in the direction from the
    /// Plane towards the Line.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// // The x axis belongs to the xy Plane
    /// let xy = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), xy.displacement(x));
    /// 
    /// // The x axis is parallel to the Plane parallel to xy
    /// let parallel_xy = Plane::new(Point::new(0.0, 0.0, -2.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(Vector::new(0.0, 0.0, 2.0), parallel_xy.displacement(x));
    /// 
    /// // The x axis intersects with the yz Plane in exactly one Point
    /// let yz = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), yz.displacement(x));
    /// ```
    fn displacement(self, line: Line) -> Vector {
        match self.relative_position(line) {
            RelativePosition::Parallel => return self.point.displacement(line),
            _ => return Vector::new(0.0, 0.0, 0.0)
        }
    }

    /// Computes the distance between a Plane and a [`Line`]. 0.0 is returned
    /// when their [`RelativePosition`] is Equivalent of Parallel.
    /// 
    /// Note: this method is a wrapper around the `Line::distance()` method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// // The x axis belongs to the xy Plane
    /// let xy = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(0.0, xy.distance(x));
    /// 
    /// // The x axis is parallel to the Plane parallel to xy
    /// let parallel_xy = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(2.0, parallel_xy.distance(x));
    /// 
    /// // The x axis intersects with the yz Plane in exactly one Point
    /// let yz = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(0.0, yz.distance(x));
    /// ```
    fn distance(self, line: Line) -> f64 {
        line.distance(self)
    }

    /// Determines the relative position of a Plane and a [`Line`].
    /// 
    /// A Line and a Plane can adopt these [`RelativePosition`]s in 3D space:
    ///  1. Equivalent
    ///  2. Parallel
    ///  3. Intersecting
    /// 
    /// A Plane and a Line are Equivalent if the Line belongs to the Plane,
    /// that is the Line's directional [`Vector`] and the Plane's normal Vector
    /// are perpendicular and the [`Point`] defining the Line lies on the Plane.
    /// 
    /// A Plane and a Line are Parallel if they do not share any Points, i.e.
    /// if their Vectors are perpendicular but the Line's Point does not lie on
    /// the Plane.
    /// 
    /// A Plane and a Line are Intersecting if they share exactly one Point. 
    /// That is the case when their Vectors are not perpendicular.
    /// 
    /// Note: this method is a wrapper around the `Line::relative_position(plane: Plane)`
    /// method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    /// 
    /// // The x axis belongs to the xy Plane
    /// let xy = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(RelativePosition::Equivalent, x.relative_position(xy));
    /// 
    /// // The x axis is parallel to the Plane parallel to xy
    /// let parallel_xy = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(RelativePosition::Parallel, x.relative_position(parallel_xy));
    /// 
    /// // The x axis intersects with the yz Plane in exactly one Point
    /// let yz = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(RelativePosition::Intersecting, x.relative_position(yz));
    /// ```
    fn relative_position(self, line: Line) -> RelativePosition {
        line.relative_position(self)
    }
}

/// Used for determining the relative position of two Planes.
impl RelativePositionTrait<Plane> for Plane {
    /// Computes the displacement between two Planes. A zero [`Vector`] is returned
    /// when their [`RelativePosition`] is Equivalent or Intersecting.
    /// 
    /// **Note:** The computed displacement Vector is in the direction from the Plane
    /// that this method is called on towards the Plane passed in to the method.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// // Identical Planes
    /// let xy = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), xy.displacement(xy));
    /// 
    /// // Parallel Planes
    /// let parallel = Plane::new(Point::new(5.0, -4.5, 9.0), Vector::new(0.0, 0.0, 5.0)).unwrap();
    /// assert_eq!(Vector::new(0.0, 0.0, 9.0), xy.displacement(parallel));
    /// 
    /// // Intersecting Planes
    /// let xz = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 3.0, 0.0)).unwrap();
    /// assert_eq!(Vector::new(0.0, 0.0, 0.0), xy.displacement(xz));
    /// ```
    fn displacement(self, other: Plane) -> Vector {
        match self.relative_position(other) {
            RelativePosition::Parallel => return self.point.displacement(other),
            _ => return Vector::new(0.0, 0.0, 0.0)
        }
    }

    /// Computes the distance between two Planes. 0.0 is returned when their
    /// [`RelativePosition`] is Equivalent or Intersecting.
    /// 
    /// Note: unlike `Plane::distance(point: Point)`, this method always returns
    /// positive distance.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// // Identical Planes
    /// let xy = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(0.0, xy.distance(xy));
    /// 
    /// // Parallel Planes
    /// let parallel = Plane::new(Point::new(5.0, -4.5, 9.0), Vector::new(0.0, 0.0, 5.0)).unwrap();
    /// assert_eq!(9.0, xy.distance(parallel));
    /// 
    /// // Intersecting Planes
    /// let xz = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 3.0, 0.0)).unwrap();
    /// assert_eq!(0.0, xy.distance(xz));
    /// ```
    fn distance(self, other: Plane) -> f64 {
        match self.relative_position(other) {
            RelativePosition::Parallel => return self.point.distance(other).abs(),
            _ => return 0.0
        }
    }

    /// Determines the relative position of two Planes.
    /// 
    /// Two Planes can adopt the following [`RelativePosition`]s in 3D space:
    ///  1. Equivalent
    ///  2. Parallel
    ///  3. Intersecting
    /// 
    /// Two Planes are Equivalent when their normal [`Vector`]s are k multiples
    /// of one another and they share their defining [`Point`]s.
    /// 
    /// They are Parallel when their normal Vectors are k multiple of one another
    /// but they do not share their defining Points.
    /// 
    /// They are Intersecting when their normal Vectors are not k multiples.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// // Identical Planes
    /// let xy = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)).unwrap();
    /// assert_eq!(RelativePosition::Equivalent, xy.relative_position(xy));
    /// 
    /// // Parallel Planes
    /// let parallel = Plane::new(Point::new(5.0, -4.5, 9.0), Vector::new(0.0, 0.0, 5.0)).unwrap();
    /// assert_eq!(RelativePosition::Parallel, xy.relative_position(parallel));
    /// 
    /// // Intersecting Planes
    /// let xz = Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 3.0, 0.0)).unwrap();
    /// assert_eq!(RelativePosition::Intersecting, xy.relative_position(xz));
    /// ```
    fn relative_position(self, other: Plane) -> RelativePosition {
        if self.normal.is_k_multiple(other.normal) {
            if self.has_point(other.point) {
                return RelativePosition::Equivalent
            } else {
                return RelativePosition::Parallel
            }
        } else {
            return RelativePosition::Intersecting
        }
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

    /// Creates a new Plane from multiple [`Point`]s, as long as a valid plane exists.
    /// Any number of Points can be passed in to this method and it will attempt to 
    /// find the best Plane to encompass them all. This is achieved by computing the
    /// centroid of the provided points, and then computing the normal [`Vector`] of each
    /// combination of 2 Points and the centroid. These normal Vectors are then averaged
    /// to obtain the Vector that is used as the normal Vector of the returned Plane.
    /// 
    /// This method can only fail if the averaged normal Vector has its components along
    /// each of the axes (x, y, z) smaller than the tolerance, in which case None is returned.
    /// 
    /// # Example
    /// 
    /// ```
    /// use point_group::geometry::*;
    /// 
    /// let p1 = Point::new(0.0, 0.0, 0.0);
    /// 
    /// let p2 = Point::new(1.0, 0.0, 0.0);
    /// let p2b = Point::new(2.0, 0.0, -0.0);
    /// 
    /// let p3 = Point::new(-1.0, 0.0, 0.0);
    /// let p3b = Point::new(-2.0, 0.0, -0.0);
    /// 
    /// let p4 = Point::new(0.0, 1.0, 0.0);
    /// let p4b = Point::new(0.0, 2.0, -0.0);
    /// 
    /// let p5 = Point::new(0.0, -1.0, 0.0);
    /// let p5b = Point::new(0.0, -2.0, -0.0);
    /// 
    /// let points = vec![&p1, &p2, &p2b, &p3, &p3b, &p4, &p4b, &p5, &p5b];
    /// 
    /// assert_eq!(Plane::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)),
    ///     Plane::from_multiple_points(&points, 1e-10));
    /// ```
    pub fn from_multiple_points(points: &Vec<&Point>, tolerance: f64) -> Option<Plane> {
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

    #[rstest]
    #[case::equivalent(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)), 
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(-5.9, 0.0, 0.0)), RelativePosition::Equivalent)]
    #[case::parallel(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)), 
        Line::new(Point::new(0.0, 0.0, 1.0), Vector::new(-5.9, 0.0, 0.0)), RelativePosition::Parallel)]
    #[case::intersecting_x_y(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, 6.0, 0.0), Vector::new(0.0, 1.0, 0.0)), RelativePosition::Intersecting)]
    #[case::intersecting_x_z(Line::new(Point::new(-1.5, 0.0, 0.0), Vector::new(0.1, 0.0, 0.0)),
        Line::new(Point::new(0.0, 0.0, 1.0), Vector::new(0.0, 0.0, 0.1)), RelativePosition::Intersecting)]
    #[case::intersecting_y_z(Line::new(Point::new(0.0, -9.9, 0.0), Vector::new(0.0, -2.2, 0.0)), 
        Line::new(Point::new(0.0, 0.0, -1.0), Vector::new(0.0, 0.0, 0.5)), RelativePosition::Intersecting)]
    #[case::intersecting_x_xy(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, 6.0, 0.0), Vector::new(1.0, 1.0, 0.0)), RelativePosition::Intersecting)]
    #[case::intersecting_x_xz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)), RelativePosition::Intersecting)]
    #[case::intersecting_x_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Intersecting)]
    #[case::intersecting_y_xy(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
        Line::new(Point::new(0.0, 6.0, 0.0), Vector::new(1.0, 1.0, 0.0)), RelativePosition::Intersecting)]
    #[case::intersecting_y_xz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)), RelativePosition::Intersecting)]
    #[case::intersecting_y_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
        Line::new(Point::new(0.0, 6.0, 0.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Intersecting)]
    #[case::intersecting_z_xy(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)),
        Line::new(Point::new(0.0, 0.0, 4.0), Vector::new(1.0, 1.0, 0.0)), RelativePosition::Intersecting)]
    #[case::intersecting_z_xz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)), RelativePosition::Intersecting)]
    #[case::intersecting_z_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)),
        Line::new(Point::new(0.0, 6.0, 0.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Intersecting)]
    #[case::intersecting_x_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 1.0, 4.0)), RelativePosition::Intersecting)]
    #[case::intersecting_y_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 1.0, 1.0)), RelativePosition::Intersecting)]
    #[case::intersecting_z_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 1.0, 4.0)), RelativePosition::Intersecting)]
    #[case::intersecting_xy_xz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 2.0, 0.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)), RelativePosition::Intersecting)]
    #[case::intersecting_xy_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, -9.0, 0.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Intersecting)]
    #[case::intersecting_xz_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 7.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Intersecting)]
    #[case::intersecting_xy_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 2.0, 0.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, -2.0, 1.0)), RelativePosition::Intersecting)]
    #[case::intersecting_xz_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 7.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(-1.0, 1.0, 5.0)), RelativePosition::Intersecting)]
    #[case::intersecting_yz_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 3.0, 7.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(-6.0, 1.0, 5.0)), RelativePosition::Intersecting)]
    #[case::intersecting_xyz_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(-8.0, 3.0, 7.0)),
        Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(-6.0, 1.0, 5.0)), RelativePosition::Intersecting)]
    #[case::skew_x_y(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, 6.0, 2.0), Vector::new(0.0, 1.0, 0.0)), RelativePosition::Skew)]
    #[case::skew_x_z(Line::new(Point::new(-1.5, 0.0, 0.0), Vector::new(0.1, 0.0, 0.0)),
        Line::new(Point::new(0.0, 4.0, 1.0), Vector::new(0.0, 0.0, 0.1)), RelativePosition::Skew)]
    #[case::skew_y_z(Line::new(Point::new(0.0, -9.9, 0.0), Vector::new(0.0, -2.2, 0.0)), 
        Line::new(Point::new(7.0, 0.0, -1.0), Vector::new(0.0, 0.0, 0.5)), RelativePosition::Skew)]
    #[case::skew_x_xy(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, 6.0, -9.0), Vector::new(1.0, 1.0, 0.0)), RelativePosition::Skew)]
    #[case::skew_x_xz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, 0.7, 0.0), Vector::new(1.0, 0.0, 1.0)), RelativePosition::Skew)]
    #[case::skew_x_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, 8.0, 0.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Skew)]
    #[case::skew_y_xy(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
        Line::new(Point::new(0.0, 6.0, -6.0), Vector::new(1.0, 1.0, 0.0)), RelativePosition::Skew)]
    #[case::skew_y_xz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
        Line::new(Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)), RelativePosition::Skew)]
    #[case::skew_y_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
        Line::new(Point::new(2.0, 6.0, 0.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Skew)]
    #[case::skew_z_xy(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)),
        Line::new(Point::new(3.0, 0.0, 0.0), Vector::new(1.0, 1.0, 0.0)), RelativePosition::Skew)]
    #[case::skew_z_xz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)),
        Line::new(Point::new(0.0, -9.0, 0.0), Vector::new(1.0, 0.0, 1.0)), RelativePosition::Skew)]
    #[case::skew_z_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)),
        Line::new(Point::new(4.0, 6.0, 0.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Skew)]
    #[case::skew_x_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
        Line::new(Point::new(0.0, -8.0, 0.0), Vector::new(1.0, 1.0, 4.0)), RelativePosition::Skew)]
    #[case::skew_y_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
        Line::new(Point::new(5.0, 0.0, 0.0), Vector::new(1.0, 1.0, 1.0)), RelativePosition::Skew)]
    #[case::skew_z_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0)),
        Line::new(Point::new(0.0, -7.0, 0.0), Vector::new(1.0, 1.0, 4.0)), RelativePosition::Skew)]
    #[case::skew_xy_xz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 2.0, 0.0)),
        Line::new(Point::new(0.0, -6.0, 0.0), Vector::new(1.0, 0.0, 1.0)), RelativePosition::Skew)]
    #[case::skew_xy_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, -9.0, 0.0)),
        Line::new(Point::new(6.0, 0.0, 0.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Skew)]
    #[case::skew_xz_yz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 7.0)),
        Line::new(Point::new(0.0, 0.0, -1.0), Vector::new(0.0, 1.0, 5.0)), RelativePosition::Skew)]
    #[case::skew_xy_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 2.0, 0.0)),
        Line::new(Point::new(7.0, 0.0, 0.0), Vector::new(1.0, -2.0, 1.0)), RelativePosition::Skew)]
    #[case::skew_xz_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 7.0)),
        Line::new(Point::new(8.0, 0.0, 0.0), Vector::new(-1.0, 1.0, 5.0)), RelativePosition::Skew)]
    #[case::skew_yz_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 3.0, 7.0)),
        Line::new(Point::new(9.0, 0.0, 0.0), Vector::new(-6.0, 1.0, 5.0)), RelativePosition::Skew)]
    #[case::skew_xyz_xyz(Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(-8.0, 3.0, 7.0)),
        Line::new(Point::new(0.0, -5.0, 0.0), Vector::new(-6.0, 1.0, 5.0)), RelativePosition::Skew)]
    fn test_line_line_relative_positions(#[case] line: Line, #[case] other: Line, #[case] expected: RelativePosition) {
        assert_eq!(expected, line.relative_position(other));
    }
}