#ifndef PHOTONS_HPP
#define PHOTONS_HPP

namespace xrp
{

    struct Photon
    {
        Photon(const xru::Point3D& origin, const xru::Vector3D& direction, const float energy, const int lerp_energy_index, const float lerp_percentage);

        float energy_;
        xru::Vector3D direction_;
        xru::Point3D origin_;

        // For lerping:
        int lerp_energy_index_;
        float lerp_percentage_;
    };

    class Source
    {
        protected:
        Source(xru::Point3D origin, const int spectrum);
        static Source* source_;

        xru::Vector3D generate_random_vector() const;
        void generate_energy(int& index, double& lerp_perc, float& energy) const;

        public:
        virtual Photon* generate_particle() const = 0;

        const std::vector<double>* spectrum_;
        xru::Point3D origin_;
    };

    class PointSource : public Source
    {
        private:
        PointSource(xru::Point3D origin, const int spectrum);

        public:
        PointSource* create_source(xru::Point3D origin, const int spectrum);
        Photon* generate_particle() const override;
    };

    class FiniteSource : public Source
    {
        private:
        FiniteSource(xru::Point3D origin, xru::Vector3D normal_, xru::Vector3D local_y, float dx, float dy, const int spectrum);

        public:
        FiniteSource* create_source(xru::Point3D origin, xru::Vector3D normal_, xru::Vector3D local_y, float dx, float dy, const int spectrum);
        Photon* generate_particle() const override;

        xru::Vector3D normal_, local_y_;
        // Half-lengths per local axis:
        float dx_, dy_;
    };

}

#endif