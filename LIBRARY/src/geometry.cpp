#include "XRaySimulator.hpp"
#include "geometry.hpp"

namespace xrg
{

    unsigned int Body::body_count = 0;

    Body::Body(const double x, const double y, const double z, int material):
        centre_(x, y, z), material_(material), absorbed_energy_(0)
    {
        body_index = body_count++;
    };

    double Body::calculate_dose() const
    {
        if (material_ == xrc::VACUUM) return 0;
        return absorbed_energy_ / (volume() * xrc::material_densities[material_]) * 1.602176565e-13; //1000(keV) / 0.001(g) * e
    }

    void xrg::Body::print(std::ostream &where) const
    {
        where << "Body " << body_index <<" at " << centre_ << ", with material " << material_;
    }

    void xrg::Body::set_material(const int material)
    {
        material_ = material;
    }

    void Body::absorb_energy(const double energy)
    {
        absorbed_energy_ += energy;
    }

    Ellipsoid::Ellipsoid(const double a, const double b, const double c, const double x, const double y, const double z) : Body(x, y, z), a_(a), b_(b), c_(c)
    {
        // This MUST be an ellipsoid
        assert((a_ != 0));
        assert((b_ != 0));
        assert((c_ != 0));
    }

    void xrg::Ellipsoid::intersect(const xrp::Photon &photon, double *intersections, int& numintersections) const
    {
        xru::QuadraticCoef* qc = intersect_coefs(photon);
        xru::QuadraticSolver(*qc, intersections, numintersections);
        delete qc;
    }

    xru::QuadraticCoef *Ellipsoid::intersect_coefs(const xrp::Photon &photon) const
    {
        __m256d mm_ph_dir, firsts, seconds, coefs, coefs2, result1, result2, subresult;
        double* result;

        mm_ph_dir = _mm256_set_pd(photon.direction_.dx, photon.direction_.dy, photon.direction_.dz, 0);
        coefs = _mm256_set_pd(a_, b_, c_, 0);
        result1 = _mm256_div_pd(mm_ph_dir, coefs);
        result2 = _mm256_mul_pd(result1, result1);
        result = (double*) &result2;

        // Assuming this is non-zero:
        auto divider = 1 / (result2[3] + result2[2] + result2[1]);

        xru::QuadraticCoef* qc = new xru::QuadraticCoef();
        firsts = _mm256_set_pd(photon.origin_.dx, photon.origin_.dy, photon.origin_.dz, 0);
        seconds = _mm256_set_pd(centre_.dx, centre_.dy, centre_.dz, 0);
        // photon.origin_.dx - centre_.dx, photon.origin_.dy - centre_.dy, photon.origin_.dz - centre_.dz
        subresult = _mm256_sub_pd(firsts, seconds);
        // (photon.origin_.dx - centre_.dx)*(photon.origin_.dx - centre_.dx),
        // (photon.origin_.dy - centre_.dy)*(photon.origin_.dy - centre_.dy),
        // (photon.origin_.dz - centre_.dz)*(photon.origin_.dz - centre_.dz)
        result1 = _mm256_mul_pd(subresult, subresult);
        // a_*a_, b_*b_, c_*c_
        coefs2 = _mm256_mul_pd(coefs, coefs);
        //(photon.origin_.dz - centre_.dz)*(photon.origin_.dz - centre_.dz) / (c_*c_),
        //(photon.origin_.dy - centre_.dy)*(photon.origin_.dy - centre_.dy) / (b_*b_),
        //(photon.origin_.dx - centre_.dx)*(photon.origin_.dx - centre_.dx) / (a_*a_)
        result2 = _mm256_div_pd(result1, coefs2);
        result = (double*) &result2;

        qc->q = ( result[3] + result[2] + result[1] - 1 )* divider;

        result1 = _mm256_mul_pd(subresult, mm_ph_dir);
        // ((photon.origin_.dz-centre_.dz)*photon.direction_.dz/(c_*c_)),
        // ((photon.origin_.dy-centre_.dy)*photon.direction_.dy/(b_*b_)),
        // ((photon.origin_.dx-centre_.dx)*photon.direction_.dx/(a_*a_))
        result2 = _mm256_div_pd(result1, coefs2);
        result = (double*) &result2;

        qc->phalf = ( result[3] + result[2] + result[1] ) * divider;

        return qc;
    }

    double Ellipsoid::volume() const
    {
        return 4/3 * a_ * b_ * c_ * xrc::PI;
    }

    void xrg::Ellipsoid::print(std::ostream &where) const
    {
        where << "Ellipsoid " << body_index << " at " << centre_ << " with half-axes " << a_ << ", " << b_ << " and " << c_
            << " and material " << material_;
    }

    Sphere::Sphere(const double r, const double x, const double y, const double z):
        Body(x, y, z), r_(r) { };

    void xrg::Sphere::intersect(const xrp::Photon &photon, double *intersections, int& numintersections) const
    {
        xru::QuadraticCoef* qc = intersect_coefs(photon);
        xru::QuadraticSolver(*qc, intersections, numintersections);
        delete qc;
    }

    xru::QuadraticCoef *Sphere::intersect_coefs(const xrp::Photon &photon) const
    {
        // Assume line's direction vector is normalized and non-zero
        auto delta = photon.origin_ - centre_;

        xru::QuadraticCoef* qc = new xru::QuadraticCoef();
        qc->q = delta.dot(delta) - r_*r_;
        qc->phalf = photon.direction_.dot(delta);

        return qc;
    }

    double Sphere::volume() const
    {
        return 4/3 * r_ * r_ * r_ * xrc::PI;
    }

    void xrg::Sphere::print(std::ostream &where) const
    {
        where << "Sphere " << body_index << " at " << centre_ << " with radius " << r_
            << " and material " << material_;
    }

    Cylinder::Cylinder(const double r, const double h, const double x, const double y, const double z): Body(x, y, z), r_(r), h_(h) { };

    void xrg::Cylinder::intersect(const xrp::Photon &photon, double *intersections, int& numintersections) const
    {
        bool intersects_planar_face = false;
        enum {UPPER = 1, LOWER = -1};

        numintersections = 0;

        // Cylinders are always along z-axis, so we need to check circle intersection, and z-bounds.
        xru::QuadraticCoef* qc = intersect_coefs(photon);
        int numroots;
        double roots[2];
        xru::QuadraticSolver(*qc, roots, numroots);

        for(int i = 0; i < numroots; i++)
            if ((photon.direction_.dz*roots[i] + photon.origin_.dz > centre_.dz - h_*0.5 - xrc::tolerance) &&
                (photon.direction_.dz*roots[i] + photon.origin_.dz < centre_.dz + h_*0.5 + xrc::tolerance))
                {
                    intersections[i] = roots[i];
                    numintersections++;
                }

        intersects_planar_face = planar_face_intersect(photon, intersections, numintersections, UPPER);
        intersects_planar_face = planar_face_intersect(photon, intersections, numintersections, LOWER) || intersects_planar_face;

        if (intersects_planar_face && numintersections == 2 && intersections[0] > intersections[1])
        {
            intersections[0] = intersections[0] + intersections[1];
            intersections[1] = intersections[0] - intersections[1];
            intersections[0] = intersections[0] - intersections[1];
        }

        delete qc;
    }

    xru::QuadraticCoef *Cylinder::intersect_coefs(const xrp::Photon &photon) const
    {
        // Cylinder is always along z-axis, so this is checking circle
        xru::Vector3D ray_projection = xru::Vector3D(photon.direction_.dx, photon.direction_.dy, 0);
        xru::QuadraticCoef* qc = new xru::QuadraticCoef();

        double divider = ray_projection.dot(ray_projection);
        divider = 1/divider;

        qc->phalf = (ray_projection.dot(photon.origin_) - ray_projection.dot(centre_)) * divider;
        qc->q = ((photon.origin_.dx - centre_.dx) * (photon.origin_.dx - centre_.dx) + 
                (photon.origin_.dy - centre_.dy) * (photon.origin_.dy - centre_.dy)
                - r_*r_) * divider;

        return qc;
    }

    bool xrg::Cylinder::planar_face_intersect(const xrp::Photon &photon, double *intersections, int &numintersections, const int side) const
    {
        xru::Vector3D face_normal = xru::Vector3D(0, 0, 1) * side;
        xru::Vector3D face_centre = xru::Vector3D(0, 0, h_*0.5) * side + centre_;

        double d = (face_centre - photon.origin_).dot(face_normal) / photon.direction_.dot(face_normal);
        if (d < 0) return false;

        if ((photon.direction_*d - face_centre).dot(photon.direction_*d - face_centre) < r_*r_)
        {
            intersections[numintersections] = d;
            numintersections++;
            return true;
        }

        return false;
    }

    double Cylinder::volume() const
    {
        return r_ * r_ * xrc::PI * h_;
    }

    void xrg::Cylinder::print(std::ostream &where) const
    {
        where << "Cylinder " << body_index << " at " << centre_ << " with radius " << r_ << " and height " << h_
            << " and material " << material_;
    }

    std::ostream& operator<<(std::ostream &out, const Body &o)
    {
        o.print(out);
        return out;
    }

}
