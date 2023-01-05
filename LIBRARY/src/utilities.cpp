#include "XRaySimulator.hpp"
#include "utilities.hpp"

namespace xru {

    Vector3D::Vector3D(): dx(0), dy(0), dz(0) { };

    Vector3D::Vector3D(const double x, const double y, const double z): dx(x), dy(y), dz(z) { };

    Vector3D::Vector3D(const Vector3D& other): dx(other.dx), dy(other.dy), dz(other.dz) { };

    Vector3D Vector3D::cross(const Vector3D &other) const
    {
        __m256d firsts = _mm256_set_pd(dy, dz, dz, dx);
        __m256d seconds = _mm256_set_pd(other.dz, other.dy, other.dx, other.dz);
        __m256d mult_result = _mm256_mul_pd(firsts, seconds);

        __m256d firsts2 = _mm256_set_pd(dx, dy, 0, 0);
        __m256d seconds2 = _mm256_set_pd(other.dy, other.dx, 0, 0);
        __m256d mult_result2 = _mm256_mul_pd(firsts2, seconds2);

        double* result = (double*) &mult_result;
        double* result2 = (double*) &mult_result2;

        // __m256d registers INVERSE ORDER OF PARAMETERS!!!
        return Vector3D(result[3] - result[2], result[1] - result[0], result2[3] - result2[2]);
    }

    double Vector3D::dot(const Vector3D &other) const
    {
        __m256d firsts = _mm256_set_pd(dx, dy, dz, 0);
        __m256d seconds = _mm256_set_pd(other.dx, other.dy, other.dz, 0);
        __m256d mult_result = _mm256_mul_pd(firsts, seconds);
        double* result = (double*) &mult_result;

        return result[3] + result[2] + result[1];
    }

    double Vector3D::euclidian_distance(const Vector3D &other) const
    {
        return abs(*this-other);
    }

    Vector3D Vector3D::operator + (Vector3D const &other) const
    {
        __m256d firsts = _mm256_set_pd(dx, dy, dz, 0);
        __m256d seconds = _mm256_set_pd(other.dx, other.dy, other.dz, 0);
        __m256d add_result = _mm256_add_pd(firsts, seconds);
        double* result = (double*) &add_result;

        return Vector3D(result[3], result[2], result[1]); 
    }

    Vector3D Vector3D::operator - (Vector3D const &other) const
    {
        __m256d firsts = _mm256_set_pd(dx, dy, dz, 0);
        __m256d seconds = _mm256_set_pd(other.dx, other.dy, other.dz, 0);
        __m256d sub_result = _mm256_sub_pd(firsts, seconds);
        double* result = (double*) &sub_result;

        return Vector3D(result[3], result[2], result[1]); 
    }

    Vector3D Vector3D::operator * (double const & scalar) const
    {
        __m256d firsts = _mm256_set_pd(dx, dy, dz, 0);
        __m256d seconds = _mm256_set_pd(scalar, scalar, scalar, 0);
        __m256d mul_result = _mm256_mul_pd(firsts, seconds);
        double* result = (double*) &mul_result;

        return Vector3D(result[3], result[2], result[1]);
    }

    Vector3D Vector3D::operator / (double const & scalar) const
    {
        return *this * (1/scalar);
    }

    Vector3D Vector3D::normed() const
    {
        double A = Vector3D::abs(*this);
        if (A <= 0) return Vector3D();
        A = 1/A;
        return *this * A;
    }

    void Vector3D::norm()
    {
        double A = Vector3D::abs(*this);
        if (A <= 0) return;
        A = 1/A;
        dx = dx * A;
        dy = dy * A;
        dz = dz * A;
    }

    double Vector3D::abs(Vector3D const &vec)
    {
        return std::sqrt(vec.dot(vec));
    }

    std::ostream& operator<<(std::ostream& os, const Vector3D& obj)
    {
        os << "(" << obj.dx << ", " << obj.dy << ", " << obj.dz << ")";
        return os;
    }


    void QuadraticSolver(QuadraticCoef const &coef, double *roots, int &numroots)
    {
        // Numercal stability possibly a problem.
        numroots     = 0;
        double delta = coef.phalf * coef.phalf - coef.q;
        if (delta < xrc::tolerance) return;

        delta           = std::sqrt(delta);
        roots[numroots] = -coef.phalf - delta;
        numroots += (roots[numroots] > xrc::tolerance);

        roots[numroots] = -coef.phalf + delta;
        numroots += (roots[numroots] > xrc::tolerance);
    }

    int index_finder(double element_value, const std::vector<double>& container)
    {
        for (size_t i=0; i < container.size(); i++)
            if (element_value <= container[i]) return i;
        return -1;
    }

    RandomGenerator* RandomGenerator::instance_ = nullptr;
    const double RandomGenerator::max_inverse = 1./RAND_MAX;

    RandomGenerator::RandomGenerator()
    {
        std::srand(std::time(nullptr));
        instance_ = this;
    }

    RandomGenerator* RandomGenerator::initialize_random_generator()
    {
        if (instance_ == nullptr)
            return new RandomGenerator();
        else return instance_;
    }

    RandomGenerator::~RandomGenerator()
    {
        RandomGenerator::instance_ = nullptr;
    }

    double xru::RandomGenerator::random_scalar()
    {
        return rand() * max_inverse * 2 - 1.;
    }

    double RandomGenerator::random_positive()
    {
        return rand() * max_inverse;
    }

    double RandomGenerator::random_range(const double min, const double max)
    {
        return rand() * max_inverse * (max - min) + min;
    }
}