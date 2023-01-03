#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

namespace xrg {

    class Body {
        public:
        Body(const double x, const double y, const double z, int material = xrc::materials::VACUUM);

        virtual void intersect(const xrp::Photon &photon, double* intersections, int& numintersections) const = 0;
        virtual xru::QuadraticCoef* intersect_coefs(const xrp::Photon &photon) const = 0;
        virtual void print(std::ostream& where) const;

        void set_material(const int material);
        void absorb_dose(const double energy);

        int material_;
        double dose_;
        unsigned int body_index;
        static unsigned int body_count;

        protected:
        xru::Point3D centre_;
    };

    std::ostream& operator<<(std::ostream& out, const Body& o);

    class Ellipsoid: public Body {
        public:
        Ellipsoid(const double a, const double b, const double c, const double x, const double y, const double z);

        virtual void intersect(const xrp::Photon &photon, double* intersections, int& numintersections) const override;
        virtual xru::QuadraticCoef* intersect_coefs(const xrp::Photon &photon) const override;

        virtual void print(std::ostream& where) const override;

        private:
        double a_, b_, c_;
    };

    class Sphere: public Body {
        public:
        Sphere(const double r, const double x, const double y, const double z);

        virtual void intersect(const xrp::Photon &photon, double* intersections, int& numintersections) const override;
        virtual xru::QuadraticCoef* intersect_coefs(const xrp::Photon &photon) const override;

        virtual void print(std::ostream& where) const override;

        private:
        double r_;
    };

    class Cylinder: public Body {
        public:
        Cylinder(const double r, const double h, const double x, const double y, const double z);

        virtual void intersect(const xrp::Photon &photon, double* intersections, int& numintersections) const override;
        virtual xru::QuadraticCoef* intersect_coefs(const xrp::Photon &photon) const override;
        bool planar_face_intersect(const xrp::Photon &photon, double* intersections, int& numintersections, const int side) const;

        virtual void print(std::ostream& where) const override;

        private:
        double r_, h_;
    };

    class Rectangle: public Body {
        public:
        Rectangle(const xru::Vector3D normal, const xru::Vector3D local_y, const double dx, const double dy, const double x, const double y, const double z);

        virtual void intersect(const xrp::Photon &photon, double* local_intersections, int& numintersections) const override;
        virtual xru::QuadraticCoef* intersect_coefs(const xrp::Photon &photon) const override;
        void get_local_coordinates(const xrp::Photon& photon, const double t, double* local_coordinates) const;

        virtual void print(std::ostream& where) const override;

        private:
        double dx_, dy_;
        xru::Vector3D normal_, local_y_;
    };

}



#endif