/// @file 
/// @brief Definition of a tool to print the cell output in VTK format

#ifndef __CELL_OUTPUT_VTK_HPP_INCLUDED
#define __CELL_OUTPUT_VTK_HPP_INCLUDED


#include <fstream>
#include <string>

#include "utils/array/array.hpp"
#include "utils/vec3/vec3.hpp"
#include "utils/mat3/mat3.hpp"


/// @brief Define the maximum number of digit in the number of steps
#define LOG_MAX_STEPS 9


class CommManager;


/// @brief Store cells data used in output writing
class CellOutput {
public:
  /// @brief Default constructor
  CellOutput() {}

  /// @brief Destructor (nothing to do)
  virtual ~CellOutput() {}

  /// @brief Constructor when you know the number of cells
  /// @param [in] n Number of cells
  // When you know it, it's the final on 0
  CellOutput(uint n) : id(n), r(n), p(n), t(n), d(n) {}

  void gather(CommManager* comm, CellOutput* sendBuffer, Array<int>& counts, Array<int>& disps);

  Array< uint > id; ///< Indices of the cells
  Array< vec3<double> > r; ///< Positions of the cells
  Array< mat3<double> > p; ///< Pressures of the cells
  Array< double > t; ///< Temperatures of the cells
  Array< double > d; ///< Densities of the cells

};



/// @brief Tool to write the cells in VTK format
struct CellWriterVTK {

  /// @brief Default constructor
  CellWriterVTK() {}
  /// @brief Destructor (nothing to do)
  virtual ~CellWriterVTK() {}

  /// @brief Constructor
  /// @param [in] buff The cells to be written
  /// @param [in] name Root name of the file where the cells will be written
  CellWriterVTK(CellOutput* buff, const std::string& name) 
    : rootName(name), extension(".vtk"), cells(buff) {
  }

  void setFilename(uint step, std::string& filename);

  void write(uint step);

  void writeHeader(std::ofstream& flux);
  void writePositions(std::ofstream& flux, const Array< uint >& id);
  void writePressures(std::ofstream& flux, const Array< uint >& id);
  void writeTemperatures(std::ofstream& flux, const Array< uint >& id);
  void writeDensities(std::ofstream& flux, const Array< uint >& id);

  std::string rootName; ///< Root name of the file where the cells will be written
  std::string extension; ///< Extension for the file where the cells will be written

  CellOutput* cells; ///< Storage for the particles
};

#endif // __CELL_OUTPUT_HPP_INCLUDED
