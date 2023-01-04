#include "XRaySimulator.hpp"
#include "scene.hpp"

namespace xrs
{
    
    Scene::Scene(xrp::Source* source, xrd::Detector* detector): source_(source), detector_(detector) { }

    Scene::~Scene()
    {
        delete source_;
        delete detector_;
        for (auto body_ptr: bodies_)
            delete body_ptr;
    }

    void Scene::add_body(xrg::Body *body)
    {
        bodies_.push_back(body);
    }


    void Scene::shoot_photons(const size_t number_of_photons)
    {
        xrp::Photon* photon;

        size_t intersection_indices[2];
        bool hit;

        std::vector<traversal_info> crossings;
        double intersections[2];
        int numintersections;

        int divider_for_progress = number_of_photons / 100;

        for (size_t i = 0; i < number_of_photons;)
        {
            photon = source_->generate_particle();
            detector_->check_hit(*photon, intersection_indices, hit);

            if (hit) [[unlikely]]
            {
                body_interactions(photon, crossings, intersections, numintersections);
                hit = check_photon_absorption(photon, crossings);
                detector_->m_[intersection_indices[0]][intersection_indices[1]]->photons += hit;
                crossings.clear();

                i++;    // Count only photons that hit detector
                if (i % divider_for_progress == 0)
                    std::cout << "\r" << int(100. * i /number_of_photons) << "\%" << std::flush;
            }

            delete photon;
        }
        std::cout << "\rDONE." << std::endl;
    }

    void Scene::body_interactions(const xrp::Photon* ph, std::vector<traversal_info>& crossings, double* intersections, int numintersections)
    {
        for (auto body: bodies_)
        {
            body->intersect(*ph, intersections, numintersections);
            if (numintersections == 2) [[unlikely]]
            {
                crossings.push_back(traversal_info(intersections[0], std::pair<xrg::Body*, bool> (body, false)));   // Entering
                crossings.push_back(traversal_info(intersections[1], std::pair<xrg::Body*, bool> (body, true)));    // Exiting
            }
        }
    }

    bool Scene::check_photon_absorption(xrp::Photon* ph, std::vector<traversal_info>& crossings)
    {
        if (crossings.empty()) return true;

        std::sort(crossings.begin(), crossings.end());
        double threshold = -log(xru::RandomGenerator::random_positive());
        
        double path_length = 0;
        double cumulative_sum = 0;
        const std::vector<double>* current_material;
        double attenuation_coef;
        std::stack<xrg::Body*> current_body;
        current_body.push(nullptr);
        std::vector<traversal_info>::iterator it1 = crossings.begin();
        std::vector<traversal_info>::iterator it2 = it1 + 1;

        for (; it2 != crossings.end(); it1++, it2++)
        {
            path_length += it2->first - it1->first;
            
            if (!it1->second.second) current_body.push(it1->second.first);
            else current_body.pop();

            if (current_body.top() == nullptr) continue;

            current_material = xrc::materials_keys.at(current_body.top()->material_);
            // ATT = ATT[i-1] + PERC * (ATT[i] - ATT[i-1])
            attenuation_coef = ph->lerp_percentage_ * current_material->at(ph->lerp_energy_index_) +
                               (1 - ph->lerp_percentage_) * current_material->at(ph->lerp_energy_index_ - 1);
            cumulative_sum += attenuation_coef * path_length;

            if (cumulative_sum > threshold)
            {
                current_body.top()->absorb_energy(ph->energy_);
                return false;
            }
        }

        return true;
    }

}
