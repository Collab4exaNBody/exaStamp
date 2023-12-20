#if __use_orchestrator

// /// @file
// /// @brief Implementation of InSituPlugin and InSituPluginArray classes


#include "orchestrator/inSituPlugins.hpp"

#include <dlfcn.h>
#include <iostream>
#include <string>
#include <dirent.h>


InSituPlugin::InSituPlugin(std::string path)
{

  // By default set as invalid until full loading
  this->is_valid = false;
  this->path = path;

  // Try to open the library
  void* new_handle = dlopen(path.c_str(), RTLD_NOW);
  
  if (!new_handle)
    {
      perror("dlopen");
      std::cout << "IP : ERROR loading plugin " << path << std::endl;
      return;
    }

  this->handle = new_handle;

  // Load the interface
  this->IP_run = (int (*)(ParticleInSitu* particles, MPI_Comm analyticsComm)) dlsym(this->handle, "IP_run");
  this->IP_name = (char * (*)())dlsym( this->handle, "IP_name" );

  if (!this->IP_run || !this->IP_name)
    {
      // The compulsory interface is not defined
      dlclose(this->handle);
      return;
    }
  else
    {
      // The interface is well defined
      this->is_valid = true;
      this->plugin_name = std::string(this->IP_name());
    }
  
}

InSituPlugin::~InSituPlugin()
{

  dlclose(this->handle);
  this->handle = NULL;
  this->path = "";
  this->plugin_id = -1;
}

InSituPluginArray::InSituPluginArray(std::string repository)
{

  this->current_idx = 0;

  // Scan the directory
  DIR * dir;
  struct dirent * file;

  dir = opendir(repository.c_str());

  if (!dir)
    {
      perror("opendir");
      std::cout << "Could not open the in situ plugin directory " << repository << std::endl;
    }

  int loaded = 0;

  while ((file = readdir(dir)) != NULL)
    {
      std::string fname(file->d_name);

      if ((fname == ".") || (fname == ".."))
        continue; // skip the directories

      if (fname.find(".so") == std::string::npos)
        continue; // skip the files that are not .so files

      // Load the plugin
      InSituPlugin* newPlugin = new InSituPlugin(repository + "/" + fname);

      if (newPlugin->valid())
        { // this is a valid in situ plugin
          newPlugin->set_id(this->current_idx++);
          this->plugins.push_back(newPlugin);
          loaded++;
        }
      else
        { // some required functions are missing
          std::cout << "The plugin " << fname << " was skipped because it was not a valid in situ plugin" << std::endl;
          delete newPlugin;
        }
    }          
}

InSituPluginArray::~InSituPluginArray()
{
  for (uint i=0; i<plugins.size(); ++i)
    {
      InSituPlugin* p = this->plugins[i];
      delete(p);
    }
  this->plugins.clear();
}

InSituPlugin* InSituPluginArray::getByName(std::string name)
{
  for (uint i=0; i<this->plugins.size(); ++i)
    {
      InSituPlugin* plug = this->plugins[i];
      if (plug->name() == name)
        return plug;
    }
  return NULL;
}

InSituPlugin* InSituPluginArray::getById(unsigned int id)
{
  for (uint i=0; i<this->plugins.size(); ++i)
    {
      InSituPlugin* plug = this->plugins[i];
      if (plug->id() == id)
        return plug;
    }
  return NULL;
}

#endif /* __use_orchestrator */
