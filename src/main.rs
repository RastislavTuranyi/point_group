#![allow(unused_variables)]
#![allow(unused_imports)]

use point_group::*;
use geometry::*;
use utils::*;
use std::f64::consts::PI;

fn main() {
    /* let l1 = Line::from_points(Point::new(2_f64.sqrt()/4.0, (2_f64.sqrt()-2.0)/4.0, 0.0), Point::new(-2_f64.sqrt()/4.0, (2.0-2_f64.sqrt())/4.0, 0.0));
    let l2 = Line::from_points(Point::new(2_f64.sqrt()/4.0, (2_f64.sqrt()-2.0)/4.0, 0.0),  Point::new(-2_f64.sqrt()/2.0, 0.0, 0.0));
    let l3 = Line::from_points(Point::new(2_f64.sqrt()/4.0, (2_f64.sqrt()-2.0)/4.0, 0.0), Point::new(-2_f64.sqrt()/4.0, (-2.0+2_f64.sqrt())/4.0, 0.0));
    
    println!("l1=l2 {} l1=l3 {}, l2=l3 {}", l1==l2, l1==l3, l2==l3); */

    /* let l1 = Line { point: Point { x: 0.0, y: 1.0, z: 0.0 }, vector: Vector { x: -3.700743415417188e-17, y: -0.2510067987608444, z: 0.0 } };
    let l2 = Line { point: Point { x: 0.0, y: 1.0, z: 0.0 }, vector: Vector { x: 0.0, y: -0.8150139559708763, z: 0.0 } };
    println!("{}, {}", l2.has_point(l1.point), l1.vector.is_k_multiple(l2.vector)); */

    /* let coordinates = vec!(Point::new(0.0, 2.0, 0.0), Point::new(2.0*(2.0 * PI / 5.0).sin(), 2.0*(2.0 * PI / 5.0).cos(), 0.0), 
    Point::new(2.0*(4.0 * PI / 5.0).sin(), -2.0*(PI / 5.0).cos(), 0.0),  Point::new(-2.0*(4.0 * PI / 5.0).sin(), -2.0*(PI / 5.0).cos(), 0.0), 
    Point::new(-2.0*(2.0 * PI / 5.0).sin(), 2.0*(2.0 * PI / 5.0).cos(), 0.0));
    //println!("{:?}", &coordinates);
    //println!("{:?}", has_plane_reflection_symmetry(&coordinates, &Plane { point: Point { x: -1.708035422500241e-17, y: -4.2700885562506023e-17, z: 0.0 }, normal: Vector { x: 0.0, y: -9.184850993605148e-17, z: -1.0 } }));
    //println!("{}", find_point_group(coordinates, 2));
    let x = Line::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));

    for point in coordinates {
        println!("{:?}", point.rotate_around_axis(x, 25.0) * 1.5)
    } */

    //println!("{}", Line { point: Point { x: -0.262190628, y: 3.78191367875, z: 2.09011912275 }, vector: Vector { x: 0.537800388, y: 0.45937334258333395, z: -0.5841187114166668 } }.approx_eq(&Line { point: Point { x: -0.262190628, y: 3.78191367875, z: 2.09011912275 }, vector: Vector { x: 0.26719927066666666, y: 0.22951216991666668, z: -0.29261673008333333 } }, 0.3));

    let coords = vec!(Point::new(0.000000000, 7.217686652, 11.428734326), Point::new(0.000000000, 7.217686652, 6.009265399), 
    Point::new(-1.121563127, 9.518432587, 9.622201030), Point::new(1.121563127, 4.916940716, 7.815798696), 
    Point::new(-1.431722865, 5.096011524, 9.622201030), Point::new(1.431722865, 9.339361779, 7.815798696), 
    Point::new(2.553285991, 7.038615844, 9.622201030), Point::new(-2.553285991, 7.396757460, 7.815798696), 
    Point::new(-0.641071760, 8.609328754, 10.884276381), Point::new(0.641071760, 5.826044549, 6.553723344), 
    Point::new(-0.884661534, 5.966681171, 10.884276381), Point::new(0.884661534, 8.468692133, 6.553723344), 
    Point::new(1.525733293, 7.077050030, 10.884276381), Point::new(-1.525733293, 7.358323273, 6.553723344), 
    Point::new(0.160642915, 9.859359996, 8.684211053), Point::new(-0.160642915, 4.576013308, 8.753788672), 
    Point::new(-2.368077682, 6.035970825, 8.684211053), Point::new(2.368077682, 8.399402479, 8.753788672), 
    Point::new(2.207434767, 5.757729134, 8.684211053), Point::new(-2.207434767, 8.677644169, 8.753788672), 
    Point::new(0.000000000, 7.217686652, 4.175528908), Point::new(1.878335430, 3.349367347, 7.231189701), 
    Point::new(2.410895073, 10.778532503, 7.231189701), Point::new(4.289230503, 6.910213198, 10.206810025), 
    Point::new(-4.289230503, 7.525160105, 7.231189701), Point::new(0.000000000, 7.217686652, 13.262470817), 
    Point::new(-1.878335430, 11.086005957, 10.206810025), Point::new(-2.410895073, 3.656840800, 10.206810025), 
    Point::new(0.919477979, 7.132157069, 3.849874402), Point::new(-0.385668198, 8.056742731, 3.849874402), 
    Point::new(-0.533809781, 6.464160155, 3.849874402), Point::new(1.221386870, 2.843768388, 6.708834522), 
    Point::new(3.177230895, 10.462397841, 6.708834522), Point::new(4.398617765, 6.088479578, 10.729165204), 
    Point::new(-4.398617765, 8.346893726, 6.708834522), Point::new(2.655922399, 3.550740987, 6.669598943), 
    Point::new(4.503629300, 7.684310088, 10.768400783), Point::new(2.163992319, 2.815619579, 8.002298073), 
    Point::new(4.894298073, 6.890725437, 9.435701652), Point::new(1.847706901, 11.351255753, 6.669598943), 
    Point::new(2.730305754, 11.292792510, 8.002298073), Point::new(-4.503629300, 6.751063216, 6.669598943), 
    Point::new(-4.894298073, 7.544647866, 8.002298073), Point::new(0.385668198, 6.378630572, 13.588125323), 
    Point::new(-0.919477979, 7.303216234, 13.588125323), Point::new(0.533809781, 7.971213149, 13.588125323), 
    Point::new(-2.655922399, 10.884632317, 10.768400783), Point::new(-2.163992319, 11.619753724, 9.435701652), 
    Point::new(-1.221386870, 11.591604915, 10.729165204), Point::new(-1.847706901, 3.084117551, 10.768400783), 
    Point::new(-2.730305754, 3.142580794, 9.435701652), Point::new(-3.177230895, 3.972975463, 10.729165204));

    //println!("{}", find_point_group(coords, 2, 0.07))
    /* let atom_1 = Point::from_array([0_f64, 0_f64, 0_f64]);
    let atom_2 = Point::from_array([1_f64, 0_f64, 0_f64]);
    let atom_3 = Point::from_array([2_f64, 0_f64, 0_f64]);
    let coordinates = vec!(atom_1, atom_2, atom_3);
    
    let result = find_point_group(coordinates, 1); */
    //print!("{:?}", result);

    //println!("{}", Line { point: Point { x: -1.3466360139166669, y: 2.7197999955000007, z: 3.592181557 }, vector: Vector { x: 0.49907284291666687, y: 0.49500360949999944, z: -0.02694136250000012 } }.approx_eq(&Line { point: Point { x: -1.3466360139166669, y: 2.7197999955000007, z: 3.592181557 }, vector: Vector { x: 1.6516119229166668, y: 1.5932588154999996, z: -0.08980453999999982 } }, 0.075));
    let coords = vec!(Point::new(5.408387870, 2.392487448, 5.455600003), Point::new(1.839182104, 1.763012487, 5.050283634), Point::new(6.793378715, 2.947163631, 5.237182883), Point::new(0.454191259, 1.208336304, 5.268700754), Point::new(6.076129674, 4.051944932, 4.343447528), Point::new(1.171440299, 0.103555003, 6.162436109), Point::new(4.816795539, 3.292984357, 4.629627883), Point::new(2.430774435, 0.862515578, 5.876255753), Point::new(3.493009179, 3.502754027, 4.174513002), Point::new(3.754560795, 0.652745908, 6.331370634), Point::new(3.006664300, 4.530076637, 3.237177884), Point::new(1.559416123, 3.833614948, 3.414097078), Point::new(5.688153851, 0.321884987, 7.091786558), Point::new(4.955958695, 1.248727753, 6.194373847), Point::new(2.291611279, 2.906772182, 4.311509790), Point::new(4.240905674, -0.374576702, 7.268705753), Point::new(7.643114730, 2.116977907, 4.594958350), Point::new(-0.395544756, 2.038522028, 5.910925287), Point::new(7.439200466, 3.390971107, 6.329479512), Point::new(-0.191630492, 0.764528829, 4.176404125), Point::new(6.459761115, 4.080950329, 3.068138157), Point::new(0.787808859, 0.074549606, 7.437745480), Point::new(6.155181702, 5.283468959, 4.817262916), Point::new(3.505241767, 4.458518844, 1.961238399), Point::new(3.098907550, 5.827174337, 3.599525698), Point::new(0.644863621, 4.621497945, 3.960718174), Point::new(1.047888354, 3.331048662, 2.289757353), Point::new(6.199681620, 0.824451273, 8.216126284), Point::new(1.092388271, -1.127969024, 5.688620720), Point::new(6.602706353, -0.465998010, 6.545165463), Point::new(4.148662424, -1.671674402, 6.906357939), Point::new(3.742328206, -0.303018909, 8.544645237));
    let cos = get_centre_of_symmetry(&coords);
    /* let planar = is_planar(&coords, 0.1);

    match planar {
        Some(plane) => println!("{:?} {:?}", plane, deviation_from_plane(&coords, plane)),
        None => println!("{:?}", deviation_from_plane(&coords, Plane::from_three_points(coords[0], coords[1], coords[3]).expect("msg")))
    } */

    //println!("{:?}", deviation_from_rotation_symmetry(&coords, Line { point: Point { x: 14.761766726, y: 15.148321191600001, z: 8.782508149 }, vector: Vector { x: 2.1966905637164165, y: 4.257967944102729, z: -2.310372869096168 } }, 4.0));

    println!("{:?}", shortest_atom_distance(&coords));

    println!("{:?}", has_planar_element(&coords, cos, 0.29));

    println!("{:?}", deviation_from_reflection_symmetry(&coords, Plane { point: Point { x: 3.6237849869375003, y: 2.0777499675000004, z: 5.252941818312499 }, normal: Vector { x: -0.17852754387578762, y: 0.5944347772145868, z: 0.7840760241937789 } }));
    //println!("{:?}", has_plane_reflection_symmetry(&coords, &Plane::new(get_centre_of_symmetry(&coords), lines[2].vector).expect("msg"), 0.2));
    //println!("{:?} {:?}", points_into_coordinates(coords), lines[2].compute_point(10.0))

    let deviations = deviation_from_plane(&coords, Plane::new(cos, Vector { x: -0.17852754387578762, y: 0.5944347772145868, z: 0.7840760241937789 }).unwrap());
    //println!("{:?} <- {:?}", deviations.iter().fold(0.0, |max, val| if val > &max {*val} else {max}), deviations);
    //println!("{:?}", deviation_from_plane(&coords, Plane::from_three_points(coords[0], coords[1], coords[3]).unwrap()));

    let vector = Plane::from_three_points(cos, coords[6], coords[8]).unwrap().normal.normalise();
    let vector = Vector { x: 0.840707187108094, y: 0.5066674461155988, z: -0.1910484875404488 };
    let deviations = deviation_from_rotation_symmetry(&coords, Line::new(cos, vector), 2.0);
    println!("{:?} <- {:?}, {:?}", deviations.iter().fold(0.0, |max, val| if val > &max {*val} else {max}), deviations, vector);
    println!("");
    /* let pp = Plane::from_three_points(cos, coords[6], coords[8]).unwrap();
    for (i, p1) in coords.iter().enumerate() {
        for p2 in coords[i+1..].iter() {
            if pp.approx_has_point(*p1, 0.3) && pp.approx_has_point(*p2, 0.3) {
                let v = Plane::from_three_points(cos, *p1, *p2).unwrap().normal.normalise();
                let devs = deviation_from_rotation_symmetry(&coords, Line::new(cos, v), 4.0);
                println!("{:?}, {:?}", devs.iter().fold(0.0, |max, val| if val > &max {*val} else {max}), v);
            }
            
        }
    } */
    //println!("{:?}", deviation_from_rotation_symmetry(&coords, Line::from_points(cos, Point::new(10.528230482, 7.339499950, 4.274270573)).unwrap(), 2.0));

    //Line { point: Point { x: 7.1446236518, y: 9.656128203933331, z: 4.7408198414000005 }, vector: Vector { x: -0.3112557392670835, y: 0.5601211141467458, z: 1.8173851621853414 } }
    /* for (i, c) in coords.iter().enumerate() {
        let v = Vector::from_two_points(cos, *c);

        /* for cc in coords[i+1..].iter() {
            let vv = Vector::from_two_points(cos, *cc);
            println!("{:?}", v.angle_degrees(vv));
        } */
        println!("{:?}", v.angle_degrees(Vector { x: 1.2527473282125494, y: -0.013485584329648914, z: -1.5910145769789892 }))
    } */
}