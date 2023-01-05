#ifndef PHOTONS_HPP
#define PHOTONS_HPP

namespace xrp
{

    struct Photon
    {
        Photon(const xru::Point3D& origin, const xru::Vector3D& direction, const float energy, const int lerp_energy_index, const float lerp_percentage);
        ~Photon() = default;

        xru::Point3D origin_;
        xru::Vector3D direction_;

        // For lerping:
        float energy_;
        int lerp_energy_index_;
        float lerp_percentage_;
    };

    std::ostream& operator<<(std::ostream& out, const Photon& p);

    class Source
    {
        protected:
        Source(const xru::Point3D& origin, const int spectrum);
        static Source* source_;

        xru::Vector3D generate_random_vector() const;
        void generate_energy(int& index, double& lerp_perc, float& energy) const;

        public:
        virtual ~Source();
        virtual Photon* generate_particle() const = 0;

        const std::vector<double>* spectrum_;
        xru::Point3D origin_;

        private:
        xru::RandomGenerator* rng_;
    };

    class PointSource : public Source
    {
        private:
        PointSource(const xru::Point3D& origin, const int spectrum);

        public:
        ~PointSource() = default;
        static PointSource* create_source(const xru::Point3D& origin, const int spectrum);
        Photon* generate_particle() const override;
    };

    class FiniteSource : public Source
    {
        private:
        FiniteSource(const xru::Point3D& origin, const xru::Vector3D& normal_, const xru::Vector3D& local_y, const float dx, const float dy, const int spectrum);

        public:
        ~FiniteSource() = default;
        static FiniteSource* create_source(const xru::Point3D& origin, const xru::Vector3D& normal_, const xru::Vector3D& local_y, const float dx, const float dy, const int spectrum);
        Photon* generate_particle() const override;

        xru::Vector3D normal_, local_y_;
        // Half-lengths per local axis:
        float dx_, dy_;
    };

}

#endif