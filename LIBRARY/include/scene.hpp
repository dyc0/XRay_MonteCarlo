#ifndef SCENE_HPP
#define SCENE_HPP

namespace xrs
{

    typedef std::pair<double, std::pair<xrg::Body*, bool>> traversal_info;

    class Scene
    {
        public:
        Scene(xrp::Source* source, xrd::Detector* detector);
        ~Scene();

        void add_body(xrg::Body* body);

        void shoot_photons(const size_t number_of_photons);
        void body_interactions(const xrp::Photon* ph, std::vector<traversal_info>& crossings, double* intersections, int numintersections);
        bool check_photon_absorption(xrp::Photon* ph, std::vector<traversal_info>& crossings, std::vector<double>& absorbed_energies);

        xrp::Source* source_;
        xrd::Detector* detector_;

        std::vector<xrg::Body*> bodies_;
        std::vector<double> absorbed_energies_;
    };

}

#endif