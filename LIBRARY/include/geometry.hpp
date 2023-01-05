#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

namespace xrg {

    class Body 
    {
        public:
        Body(const double x, const double y, const double z, int material = xrc::materials::VACUUM);
        virtual ~Body() = default;

        virtual void intersect(const xrp::Photon &photon, double* intersections, int& numintersections) const = 0;
        virtual xru::QuadraticCoef* intersect_coefs(const xrp::Photon &photon) const = 0;

        virtual double volume() const = 0;
        double calculate_dose() const;

        virtual void print(std::ostream& where) const;

        void set_material(const int material);
        void absorb_energy(const double energy);

        xru::Point3D centre_;

        int material_;
        double absorbed_energy_;
        unsigned int body_index;
        static unsigned int body_count;
    };

    std::ostream& operator<<(std::ostream& out, const Body& o);

    class Ellipsoid: public Body 
    {
        public:
        Ellipsoid(const double a, const double b, const double c, const double x, const double y, const double z);
        ~Ellipsoid() = default;

        virtual void intersect(const xrp::Photon &photon, double* intersections, int& numintersections) const override;
        virtual xru::QuadraticCoef* intersect_coefs(const xrp::Photon &photon) const override;

        virtual double volume() const override;

        virtual void print(std::ostream& where) const override;

        private:
        double a_, b_, c_;
    };

    class Sphere: public Body 
    {
        public:
        Sphere(const double r, const double x, const double y, const double z);
        ~Sphere() = default;

        virtual void intersect(const xrp::Photon &photon, double* intersections, int& numintersections) const override;
        virtual xru::QuadraticCoef* intersect_coefs(const xrp::Photon &photon) const override;

        virtual double volume() const override;

        virtual void print(std::ostream& where) const override;

        private:
        double r_;
    };

    class Cylinder: public Body 
    {
        public:
        Cylinder(const double r, const double h, const double x, const double y, const double z);
        ~Cylinder() = default;

        virtual void intersect(const xrp::Photon &photon, double* intersections, int& numintersections) const override;
        virtual xru::QuadraticCoef* intersect_coefs(const xrp::Photon &photon) const override;
        bool planar_face_intersect(const xrp::Photon &photon, double* intersections, int& numintersections, const int side) const;

        virtual double volume() const override;

        virtual void print(std::ostream& where) const override;

        private:
        double r_, h_;
    };

}



#endif