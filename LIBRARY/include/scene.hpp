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
        bool check_photon_absorption(xrp::Photon* ph, std::vector<traversal_info>& crossings);

        void save_absorptions(const std::string& path) const;

        xrp::Source* source_;
        xrd::Detector* detector_;
        std::vector<xrg::Body*> bodies_;
    };

}

#endif