// Ourselves:
#include <fitter_module.h>

// Standard library:
#include <iostream>
#include <stdexcept>

// geomtools
#include "bayeux/geomtools/line_3d.h"
#include "bayeux/geomtools/helix_3d.h"

// falaise
#include <falaise/snemo/datamodels/data_model.h>
#include <falaise/snemo/datamodels/tracker_trajectory.h>
#include <falaise/snemo/datamodels/tracker_trajectory_solution.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include "falaise/snemo/datamodels/calibrated_calorimeter_hit.h"
#include "falaise/snemo/datamodels/base_trajectory_pattern.h"

// Registration instantiation macro :
DPP_MODULE_REGISTRATION_IMPLEMENT(fitter_module,
				  "fitter_module")


void fitter_module::initialize(const datatools::properties  & setup_,
					       datatools::service_manager   & service_manager_,
					       dpp::module_handle_dict_type & /* module_dict_ */)
{
  DT_THROW_IF (this->is_initialized(),
	       std::logic_error,
	       "Module 'Fitter' is already initialized ! ");
  
  dpp::base_module::_common_initialize(setup_);
 
  // Look for services
  if (service_manager_.has("geometry"))
  {
    const geomtools::geometry_service& GS = service_manager_.get<geomtools::geometry_service> ("geometry");

    // initialize geometry manager
    //    std::cout << "Initialize geo manager " << std::endl;
    geometry_manager_ = &GS.get_geom_manager();
    DT_THROW_IF(!geometry_manager_,
                std::runtime_error,
                "Null pointer to geometry manager return by geometry_service");
  }

  // check the label
  _TTD_label_ = snemo::datamodel::data_info::default_tracker_trajectory_data_label();

  eventCounter = 0;
  this->_set_initialized(true);

  return;
}

void fitter_module::reset()
{
  DT_THROW_IF (! this->is_initialized(),
	       std::logic_error,
	       "Module 'Fitter' is not initialized !");
  this->_set_initialized(false);
  
  // clean up
  _TTD_label_.clear();

  eventCounter = 0;
  std::cout << "Fitter Module finished." << std::endl;
  return;
}

// Constructor :
fitter_module::fitter_module(datatools::logger::priority logging_priority_)
  : dpp::base_module(logging_priority_)
{
}

// Destructor :
fitter_module::~fitter_module()
{
  // MUST reset module at destruction
  if (this->is_initialized()) reset();
}

// Processing :
dpp::base_module::process_status fitter_module::process(datatools::things & data_record_)
{
  DT_THROW_IF (! this->is_initialized(), std::logic_error,
	       "Module 'Fitter' is not initialized !");
  
  ///////////////////////////////////
  // Check tracker clustering data //
  ///////////////////////////////////
  
  bool preserve_former_output = false;
  
  // check if some 'tracker_clustering_data' are available in the data model:
  snemo::datamodel::tracker_clustering_data * ptr_cluster_data = 0;
  if (! data_record_.has(_TCD_label_)) {
    ptr_cluster_data = &(data_record_.add<snemo::datamodel::tracker_clustering_data>(_TCD_label_));
  } else {
    ptr_cluster_data = &(data_record_.grab<snemo::datamodel::tracker_clustering_data>(_TCD_label_));
  }
  snemo::datamodel::tracker_clustering_data & the_clustering_data = *ptr_cluster_data;
  if (the_clustering_data.has_solutions()) 
    if (! preserve_former_output) 
      the_clustering_data.reset();


  ///////////////////////////////////
  // Check tracker trajectory data //
  ///////////////////////////////////

  //  bool preserve_former_output = true; // keep all
  
  // check if some 'tracker_trajectory_data' are available in the data model:
  // snemo::datamodel::tracker_trajectory_data * ptr_trajectory_data = 0;
  // if (! data_record_.has(_TTD_label_)) {
  //   ptr_trajectory_data = &(data_record_.add<snemo::datamodel::tracker_trajectory_data>(_TTD_label_));
  // } else {
  //   ptr_trajectory_data = &(data_record_.grab<snemo::datamodel::tracker_trajectory_data>(_TTD_label_));
  // }
  // snemo::datamodel::tracker_trajectory_data & the_trajectory_data = *ptr_trajectory_data;
  // if (the_trajectory_data.has_solutions()) 
  //   if (! preserve_former_output) 
  //     the_trajectory_data.reset();
  
  
  /********************
   * Process the data *
   ********************/
  
  // Main processing method :
  // Process the fitter :
  namespace sdm = snemo::datamodel;

  // Process events for trajectory consolidation
  std::cout << "In process: event counter = " << eventCounter << std::endl;

  // get all trajectory solutions
  // const sdm::tracker_trajectory_data::solution_col_type& all_solutions = the_trajectory_data.get_solutions();
  // for (auto entry : all_solutions) { 
  //   //entry.get().tree_dump();
  //   const sdm::tracker_trajectory_solution::trajectory_col_type &defaults = entry.get().get_trajectories();
    
  //   for (auto traj_handle : defaults) {
  //     const sdm::base_trajectory_pattern & the_base_pattern = traj_handle.get().get_pattern();
  //     if (the_base_pattern.get_pattern_id()=="line") {
  // 	const geomtools::line_3d & the_shape = (const geomtools::line_3d&)the_base_pattern.get_shape();
  // 	geomtools::vector_3d dir = the_shape.get_direction_on_curve(the_shape.get_first());
  // 	std::cout << "directional: " << dir.x() << ", " << dir.y() << ", " << dir.z() << std::endl;
  // 	std::cout << "p on line: " << the_shape.get_first().x() << ", " << the_shape.get_first().y() << ", " << the_shape.get_first().z() << std::endl;
  //     }
  //     else {
  // 	const geomtools::helix_3d & the_shape = (const geomtools::helix_3d&)the_base_pattern.get_shape();
  // 	the_shape.tree_dump();
  //     }
  //   }
  // }
  

  eventCounter++;
  return dpp::base_module::PROCESS_SUCCESS;
}

