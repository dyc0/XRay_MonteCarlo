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
        body->set_index(bodies_.size());
        bodies_.push_back(body);
        absorbed_energies_.push_back(0);
    }


    void Scene::shoot_photons(const size_t number_of_photons)
    {
        xrp::Photon* photon;

        size_t intersection_indices[2] = {0, 0};
        bool hit = false;

        std::vector<traversal_info> crossings;
        double intersections[2] = {0, 0};
        int numintersections = 0;

        int divider_for_progress;
        double progress = 0;
        int private_progress;
        std::vector<double> private_energies;
        std::vector<std::vector<size_t>> private_pixels;

        omp_set_num_threads(12);

        #pragma omp parallel private(photon, intersection_indices, hit, crossings, intersections, numintersections, private_energies, private_pixels, private_progress) shared(divider_for_progress, progress)
        {
            #pragma omp single
            {
                divider_for_progress =  number_of_photons / 100 / omp_get_num_threads();
                divider_for_progress += (divider_for_progress == 0)*100;
            }

            for (size_t i = 0; i < detector_->npx_y_; i++)
            {
                private_pixels.push_back(std::vector<size_t>());
                for (size_t j = 0; j < detector_->npx_x_; j++)
                    private_pixels.back().push_back(0);
            }

            for (auto body: bodies_)
                private_energies.push_back(0);

            private_progress = 0;

            #pragma omp for 
            for (size_t i = 0; i < number_of_photons; i++)
            {
                photon = new xrp::Photon();
                hit = false;
                while(!hit)
                {
                    delete photon;
                    photon = source_->generate_particle();
                    detector_->check_hit(*photon, intersection_indices, hit);
                }


                body_interactions(photon, crossings, intersections, numintersections);
                hit = check_photon_absorption(photon, crossings, private_energies);
                private_pixels[intersection_indices[0]][intersection_indices[1]] += hit;
                crossings.clear();
                
                private_progress++;    // Count only photons that hit detector
                if (private_progress % divider_for_progress == 0)
                    #pragma omp critical
                    {
                        progress += private_progress;
                        private_progress = 0;
                        std::cout << "\r" << int(100 * progress / number_of_photons) << "\%" << std::flush;
                    }

                delete photon;
            }

            #pragma omp single
            std::cout << "\rDONE." << std::endl;

            #pragma omp critical
            {
                std::cout << std::endl << "ENERGIES " << omp_get_thread_num() << ":";
                for (size_t i = 0; i < bodies_.size(); i++)
                {
                    bodies_[i]->absorb_energy(private_energies[i]);
                    std::cout<<private_energies[i] << " ";
                }
                for (size_t i = 0; i < detector_->npx_y_; i++)
                    for (size_t j = 0; j < detector_->npx_x_; j++)
                        detector_->m_[i][j]->photons += private_pixels[i][j];
            }
        }
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

    bool Scene::check_photon_absorption(xrp::Photon* ph, std::vector<traversal_info>& crossings, std::vector<double>& absorbed_energies)
    {
        if (crossings.empty()) return true;

        std::sort(crossings.begin(), crossings.end());
        double treshold = -log(xru::RandomGenerator::random_positive());
        
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

            if (cumulative_sum > treshold)
            {
                absorbed_energies[current_body.top()->body_index] += ph->energy_;
                return false;
            }
        }

        // Lastly, is it absorbed in detector coating?
        return detector_->coating_hit(ph, cumulative_sum, treshold);
    }

}
