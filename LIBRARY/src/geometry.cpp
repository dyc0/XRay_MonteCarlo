#include "XRaySimulator.hpp"
#include "geometry.hpp"

namespace xrg
{

    unsigned int Body::body_count = 0;

    Body::Body(const double x, const double y, const double z, int material): centre_(x, y, z), material_(material), dose_(0)
    {
        body_index = body_count++;
    };

    xru::QuadraticCoef *xrg::Body::intersect_coefs(const xrp::Photon &photon) const
    {
        return nullptr;
    }

    void xrg::Body::print(std::ostream &where) const
    {
        where << "Body " << body_index <<" at " << centre_ << ", with material " << material_;
    }

    void xrg::Body::set_material(const int material)
    {
        material_ = material;
    }

    void Body::absorb_dose(const double energy)
    {
        return;
    }

    Ellipsoid::Ellipsoid(const double a, const double b, const double c, const double x, const double y, const double z) : a_(a), b_(b), c_(c), Body(x, y, z)
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
    }

    xru::QuadraticCoef *Ellipsoid::intersect_coefs(const xrp::Photon &photon) const
    {
        // Assuming this is non-zero:
        auto divider = 1 / ((photon.direction_.dz/c_)*(photon.direction_.dz/c_) +
                            (photon.direction_.dy/b_)*(photon.direction_.dy/b_) +
                            (photon.direction_.dx/a_)*(photon.direction_.dx/a_));

        xru::QuadraticCoef* qc = new xru::QuadraticCoef();
        qc->q = ( (photon.origin_.dz - centre_.dz)*(photon.origin_.dz - centre_.dz) / (c_*c_) +
                (photon.origin_.dy - centre_.dy)*(photon.origin_.dy - centre_.dy) / (b_*b_) +
                (photon.origin_.dx - centre_.dx)*(photon.origin_.dx - centre_.dx) / (a_*a_) -
                1 )* divider;
        qc->phalf = ( ((photon.origin_.dz-centre_.dz)*photon.direction_.dz/(c_*c_)) +
                    ((photon.origin_.dy-centre_.dy)*photon.direction_.dy/(b_*b_)) +
                    ((photon.origin_.dx-centre_.dx)*photon.direction_.dx/(a_*a_)) ) * divider;

        return qc;
    }

    void xrg::Ellipsoid::print(std::ostream &where) const
    {
        where << "Ellipsoid " << body_index << " at " << centre_ << " with half-axes " << a_ << ", " << b_ << " and " << c_
            << " and material " << material_;
    }

    Sphere::Sphere(const double r, const double x, const double y, const double z):
        r_(r), Body(x, y, z) { };

    void xrg::Sphere::intersect(const xrp::Photon &photon, double *intersections, int& numintersections) const
    {
        xru::QuadraticCoef* qc = intersect_coefs(photon);
        xru::QuadraticSolver(*qc, intersections, numintersections);
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

    void xrg::Sphere::print(std::ostream &where) const
    {
        where << "Sphere " << body_index << " at " << centre_ << " with radius " << r_
            << " and material " << material_;
    }

    Cylinder::Cylinder(const double r, const double h, const double x, const double y, const double z): Body(x, y, z), h_(h), r_(r) { };

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

    Rectangle::Rectangle(const xru::Vector3D normal, const xru::Vector3D local_y, const double dx, const double dy, const double x, const double y, const double z):
        Body(x, y, z), dx_(dx), dy_(dy), normal_(normal), local_y_(local_y)
        {
            normal_.norm();
            local_y_.norm();
        }

    void Rectangle::intersect(const xrp::Photon &photon, double *local_intersections, int &numintersections) const
    {
        // CHECK: This code might be mightly inefficient.

        numintersections = 0;

        // First check plane intersection - the normal should be directed to origin
        if (photon.direction_.dot(centre_) < xrc::tolerance) [[likely]]     // -> likely because this is used in photon generation
            return;
        // Check if it's parallel
        if (abs(photon.direction_.dot(normal_)) < xrc::tolerance) [[unlikely]]
            return;

        double t = (centre_ - photon.origin_).dot(normal_) / photon.direction_.dot(normal_);
        // Check if it's behind
        if (t < xrc::tolerance) return;

        // Get local coordinates
        xru::Vector3D local_vec = photon.origin_ + photon.direction_ * t - centre_;
        local_intersections[0] = local_vec.dot(local_y_.cross(normal_));
        local_intersections[1] = local_vec.dot(local_y_);

        numintersections += (  local_intersections[0] < dx_ - xrc::tolerance && local_intersections[0] > -dx_ + xrc::tolerance
                            && local_intersections[1] < dy_ - xrc::tolerance && local_intersections[1] > -dy_ + xrc::tolerance );
    }

    void Rectangle::print(std::ostream & where) const
    {
        where << "Rectangle " << body_index << " at " << centre_ << " with half-dimensions (" << dx_ << ", " << dy_ << ") "
              << " facing to " << normal_ << " with local y-axis " << local_y_ << " and material " << material_;
    }

    xru::QuadraticCoef *Rectangle::intersect_coefs(const xrp::Photon &photon) const
    {
        return nullptr;
    }

}
