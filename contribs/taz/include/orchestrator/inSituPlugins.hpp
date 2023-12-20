#if __use_orchestrator

/// @file
/// @brief Definition of classes to hold informations about the plugins
///
///

#ifndef __INSITUPLUGINS_HPP_INCLUDED
#define __INSITUPLUGINS_HPP_INCLUDED


#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>

#include "io/particleInSitu.hpp"


/// @brief Class that defines an plugin for in situ processing
class InSituPlugin
{
  
public:

  /// @brief Default constructor
  /// @param[in] path The path were is located the .so file that defines the plugin
  InSituPlugin(std::string path);

  ~InSituPlugin();

  /// @brief Set the plugin id
  void set_id(unsigned int id) { this->plugin_id = id; }
  /// @brief Accessor to the plugin id
  unsigned int id() { return this->plugin_id; }
  /// @brief Accessor to the flag valid
  bool valid() { return this->is_valid; }
  /// @brief Accessor to the plugin name
  std::string name() { return this->plugin_name; }

  // Plugin interface
  /// @brief Function that runs the plugin
  /// @param[in] particles Pointer to the copied data
  /// @param[in] analyticsComm MPI communicator for analytics execution
  int (*IP_run) (ParticleInSitu* particles, MPI_Comm analyticsComm);

  /// @brief Function that returns the name of the plugin
  char * (*IP_name)();

private:

  std::string path;         ///< Path to the .so
  void* handle;             ///< Dlopen handle
  bool is_valid;            ///< Tells whether it is a valid plugin for in situ processing
  unsigned int plugin_id;   ///< Id for the plugin
  std::string plugin_name;  ///< Name of the plugin
  
};

/// @brief Class that defines an array of in situ plugin
class InSituPluginArray
{

public:

  InSituPluginArray();
  InSituPluginArray(std::string repository);
  ~InSituPluginArray();

  /// @brief Get a pointer to a plugin given its name
  InSituPlugin* getByName(std::string name);
  /// @brief Get a pointer to a plugin given its id
  InSituPlugin* getById(unsigned int id);

private:

  std::vector<InSituPlugin*> plugins; ///< Vector that stores pointers to the plugins
  unsigned int current_idx;           ///< Variable that stores the current number of plugins in the vector
  
};

#endif /* __INSITUPLUGINS_HPP_INCLUDED */

#endif /* __use_orchestrator */
