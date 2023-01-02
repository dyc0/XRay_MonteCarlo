#include "XRaySimulator.hpp"

using namespace xru;

int main(int argc, const char* argv[])
{   
    const double tolerance = 0.0000001;
    std::cout << "\n\t-------\n\tTESTING\n\t-------\n" << std::endl;

    Vector3D a;
    Vector3D b(1, 0, 0);
    Vector3D c(0, 1, 0);
    Vector3D d(5, 3, 1);
    Vector3D b2(b);

    std::cout << "a  = " << a  << std::endl;
    std::cout << "b  = " << b  << std::endl;
    std::cout << "c  = " << c  << std::endl;
    std::cout << "d  = " << d  << std::endl;
    std::cout << "b2 = " << b2 << std::endl;
    std::cout << std::endl << std::endl;

    std::cout << "\t----------------\n\tCHECKING VECTORS\n\t----------------\n";

    // RETURNS 1
    std::cout << "\nCHECKING BASIC VECTOR FUNCTIONS\n\n";
    std::cout << "b cross c = " << b.cross(c) << std::endl;
    if (b.cross(c).dz != 1) return 1;
    std::cout << "b dot   c = " << b.dot(c) << std::endl;
    if (b.dot(c) != 0) return 1;
    std::cout << "b Euclidian distance from c = " << b.euclidian_distance(c) << std::endl;
    if (abs(b.euclidian_distance(c) - sqrt(2)) > tolerance) return 1;
    std::cout << "b cross d = " << b.cross(d) << std::endl;
    std::cout << "b dot   d = " << b.dot(d) << std::endl;
    std::cout << "b cross a = " << b.cross(a) << std::endl;
    std::cout << "b dot   a = " << b.dot(a) << std::endl;
    std::cout << std::endl;

    // RETURNS 2
    std::cout << "\nCHECKING BASIC VECTOR OPERATORS\n\n";
    std::cout << "abs d = " << Vector3D::abs(d) << std::endl;
    if (Vector3D::abs(b) != 1) return 2;
    std::cout << "d - b = " << d-b << std::endl;
    if (Vector3D::abs(b-b2) != 0) return 2;

    // RETURNS 3
    std::cout << "\nCHECKING BASIC VECTOR NORMING\n\n";
    std::cout << "norm d = " << d.normed() << std::endl;
    std::cout << "norm c = " << c.normed() << std::endl;
    std::cout << "norm a = " << a.normed() << std::endl;
    if (abs(Vector3D::abs(d.normed()) - 1) > tolerance) return 2;
    if (Vector3D::abs(a.normed()) != 0) return 2;
    std::cout << "d = " << d << std::endl;
    d.norm();
    std::cout << "normed d = " << d << std::endl;
    if (abs(Vector3D::abs(d) - 1) > tolerance) return 2;
    std::cout << std::endl;

    // RETURNS 4
    std::cout << "\n\t-------------------------\n\tCHECKING RANDOM GENERATOR\n\t-------------------------\n\n";
    
    std::cout << "Initializing RNG." << std::endl << std::endl;
    auto rng = RandomGenerator::initialize_random_generator();
    
    double rnd_num;
    int num_cols = 10;
    
    std::cout << "Generating 100 random scalars:" << std::endl;
    for (int i = 1; i<101; i++)
    {
        rnd_num = RandomGenerator::random_scalar();
        if (rnd_num < -1 || rnd_num > 1) return 4;
        std::cout << rnd_num << " ";
        if (i % num_cols == 0) std::cout << std::endl;
    }
    std::cout<<std::endl;
    
    std::cout << "Generating 100 random positives:" << std::endl;
    for (int i = 1; i<101; i++)
    {
        rnd_num = RandomGenerator::random_positive();
        if (rnd_num < 0 || rnd_num > 1) return 4;
        std::cout << rnd_num << " ";
        if (i % num_cols == 0) std::cout << std::endl;
    }
    std::cout<<std::endl;
    
    std::cout << "Generating 100 randoms in range (-5,5):" << std::endl;
    for (int i = 1; i<101; i++)
    {
        rnd_num = RandomGenerator::random_range(-5, 5);
        if (rnd_num < -5 || rnd_num > 5) return 4;
        std::cout << rnd_num << " ";
        if (i % num_cols == 0) std::cout << std::endl;
    }
    std::cout<<std::endl;
    
    std::cout << "Generating 100 randoms in range (10,20):" << std::endl;
    for (int i = 1; i<101; i++)
    {
        rnd_num = RandomGenerator::random_range(10, 20);
        if (rnd_num < 10 || rnd_num > 20) return 4;
        std::cout << rnd_num << " ";
        if (i % num_cols == 0) std::cout << std::endl;
    }
    std::cout<<std::endl;

    delete rng;

    // RETURNS 5
    std::cout << "\n\t----------------------\n\tCHECKING INDEX FINDING\n\t----------------------\n\n";
    int index;
    index = index_finder(52.3, xrc::energies);
    std::cout << "Index for energy value 52.3 is " << index << "." << std::endl;
    if (index != 93) return 5;
    index = index_finder(6.1, xrc::energies);
    std::cout << "Index for energy value 6.1 is " << index << "." << std::endl;
    if (index != 1) return 5;
    index = index_finder(4.1, xrc::energies);
    std::cout << "Index for energy value 4.1 is " << index << "." << std::endl;
    if (index != 0) return 5;
    index = index_finder(500, xrc::energies);
    std::cout << "Index for energy value 500 is " << index << "." << std::endl;
    if (index != -1) return 5;
    std::cout << std::endl;

    return 0;
}