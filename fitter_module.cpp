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
#include <falaise/snemo/datamodels/calibrated_data.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>

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
  _TCD_label_ = snemo::datamodel::data_info::default_tracker_clustering_data_label();

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
  
  snemo::datamodel::tracker_clustering_data * ptr_cluster_data = 0;
  if (!data_record_.has("TCD")) {
    std::cerr << "failed to grab TCD bank " << std::endl;
    return dpp::base_module::PROCESS_INVALID;
  }
  else {
    ptr_cluster_data = &(data_record_.grab<snemo::datamodel::tracker_clustering_data>(_TCD_label_));
    snemo::datamodel::tracker_clustering_data & the_clustering_data = *ptr_cluster_data;
    if (! the_clustering_data.has_solutions()) {
      std::cerr << "TCD bank empty, no clustering solution, no fitting." << std::endl;
      return dpp::base_module::PROCESS_INVALID;
    }
  }


  ///////////////////////////////////
  // Check tracker trajectory data //
  ///////////////////////////////////

  bool preserve_former_output = true; // keep all
  
  // check if some 'tracker_trajectory_data' are available in the data model:
  snemo::datamodel::tracker_trajectory_data * ptr_trajectory_data = 0;
  if (! data_record_.has(_TTD_label_)) {
    ptr_trajectory_data = &(data_record_.add<snemo::datamodel::tracker_trajectory_data>(_TTD_label_));
  } else {
    ptr_trajectory_data = &(data_record_.grab<snemo::datamodel::tracker_trajectory_data>(_TTD_label_));
  }
  snemo::datamodel::tracker_trajectory_data & the_trajectory_data = *ptr_trajectory_data;
  if (the_trajectory_data.has_solutions()) 
    if (! preserve_former_output) 
      the_trajectory_data.reset();
  
  
  /********************
   * Process the data *
   ********************/
  
  // Main processing method :
  // Process the fitter :
  namespace sdm = snemo::datamodel;
  GeigerRing ring;
  MetaInfo mi;
  TrackerHit th;
  std::vector<TrackerHit> rings;
  SNFitter snf;

  // Process events for trajectory consolidation
  // make a trajectory solution
  sdm::tracker_trajectory_solution::handle_type htts(new sdm::tracker_trajectory_solution);
  the_trajectory_data.add_solution(htts, true);
  the_trajectory_data.grab_solutions().back().grab().set_solution_id(the_trajectory_data.get_number_of_solutions() - 1);
  sdm::tracker_trajectory_solution & trajectory_solution = the_trajectory_data.grab_solutions().back().grab(); // maybe store in here a bit

  // Process clusters hits for fitting
  std::cout << "In process: event counter = " << eventCounter << std::endl;

  // get all cluster solutions
  const sdm::tracker_clustering_data::solution_col_type& all_solutions = ptr_cluster_data->get_solutions();

  for (auto entry : all_solutions) { 
    const sdm::tracker_clustering_solution::cluster_col_type &defaults = entry.get().get_clusters();

    for (auto cl_handle : defaults) {
      const sdm::calibrated_tracker_hit::collection_type & gg_hits_col = cl_handle.get().get_hits();

      for (auto hit_handle : gg_hits_col) {
	// work with geiger hits as members of a given cluster
	const sdm::calibrated_tracker_hit & hit = hit_handle.get();
	ring.rerr   = hit.get_sigma_r();
	ring.zerr   = hit.get_sigma_z();
	ring.radius = hit.get_r();
	ring.wirex  = hit.get_x();
	ring.wirey  = hit.get_y();
	ring.zcoord = hit.get_z();
	mi.hitid  = hit.get_id();
	mi.side   = hit.get_side();
	mi.row    = hit.get_row();
	mi.column = hit.get_layer();
	th.mi = mi;
	th.gr = ring;
	rings.push_back(th);
      }
      // ready to fit
      snf.setData(rings);
      std::vector<HelixFit>      hres = snf.fithelix();
      std::vector<LineFit>       lres = snf.fitline();
      std::vector<BrokenLineFit> bres = snf.fitbrokenline();

      // try to do something with the results
      // Line first
      for (LineFit entry : lres) {
	if (entry.status>1)
	  std::cout << "This line fit failed with status " << entry.status  << std::endl;
	else {
	  std::cout << "Line fit: (status, chi2) " << entry.status << ", " << entry.chi2 << std::endl;
	  std::cout << "slope, intercept in xy " << entry.slxy << ", " << entry.ixy << std::endl;
	  std::cout << "errors in xy " << entry.errslxy << ", " << entry.errixy << std::endl;
	  std::cout << "slope, intercept in xz " << entry.slxz << ", " << entry.ixz << std::endl;
	  std::cout << "errors in xz " << entry.errslxz << ", " << entry.errixz << std::endl;
	}
      }

      // Helix
      for (HelixFit entry : hres) {
	if (entry.status>1)
	  std::cout << "This helix fit failed with status " << entry.status  << std::endl;
	else {
	  std::cout << "Helix fit: (status, chi2) " << entry.status << ", " << entry.chi2 << std::endl;
	  std::cout << "radius and pitch " << entry.radius << ", " << entry.pitch << std::endl;
	  std::cout << "errors in r, p " << entry.raderr << ", " << entry.errpitch << std::endl;
	  std::cout << "centre (x,y,z) " << entry.xc << ", " << entry.yc << ", " << entry.zc << std::endl;
	  std::cout << "centre errors " << entry.errxc << ", " << entry.erryc << ", " << entry.errzc << std::endl;
	}
      }

      // Broken Line last
      for (BrokenLineFit entry : bres) {
	if (entry.status>1)
	  std::cout << "This broken line fit failed with status " << entry.status  << std::endl;
	else {
	  LineFit lf = entry.linefit1;
	  if (lf.chi2<0.0) { // just one linefit, no break
	    std::cout << "Broken Line fit: (status, chi2) " << entry.status << ", " << entry.chi2 << std::endl;
	    std::cout << "slope, intercept in xy " << lf.slxy << ", " << lf.ixy << std::endl;
	    std::cout << "errors in xy " << lf.errslxy << ", " << lf.errixy << std::endl;
	    std::cout << "slope, intercept in xz " << lf.slxz << ", " << lf.ixz << std::endl;
	    std::cout << "errors in xz " << lf.errslxz << ", " << lf.errixz << std::endl;
	  }
	  else { // two lines fitted from break
	    LineFit lf2 = entry.linefit2;
	    std::cout << "Broken Line fit: (status, chi2) " << entry.status << ", " << entry.chi2 << std::endl;
	    std::cout << "(1) slope, intercept in xy " << lf.slxy << ", " << lf.ixy << std::endl;
	    std::cout << "(1) errors in xy " << lf.errslxy << ", " << lf.errixy << std::endl;
	    std::cout << "(1) slope, intercept in xz " << lf.slxz << ", " << lf.ixz << std::endl;
	    std::cout << "(1) errors in xz " << lf.errslxz << ", " << lf.errixz << std::endl;
	    std::cout << "Second line " << std::endl;
	    std::cout << "(2) slope, intercept in xy " << lf2.slxy << ", " << lf2.ixy << std::endl;
	    std::cout << "(2) errors in xy " << lf2.errslxy << ", " << lf2.errixy << std::endl;
	    std::cout << "(2) slope, intercept in xz " << lf2.slxz << ", " << lf2.ixz << std::endl;
	    std::cout << "(2) errors in xz " << lf2.errslxz << ", " << lf2.errixz << std::endl;
	  }
	}
      }
      rings.clear();
    }
  }

  eventCounter++;
  return dpp::base_module::PROCESS_SUCCESS;
}

