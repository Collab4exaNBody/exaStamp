#include "utils/array/extArray.hpp"
#include <fstream>

/// @brief A class to measure time
///
/// Measure and store time used for the different program parts of the program for each step :
class metricsDetails{

  private :

  ExtArray<double> m_time_generic;    ///< Store each time used for generic timers
  ExtArray<double> m_time_neighbours; ///< Store each time used for neighbours timers
  ExtArray<double> m_time_potential;  ///< Store each time used for potential timers
  ExtArray<double> m_time_ghost;      ///< Store each time used for ghost timers
  ExtArray<double> m_time_refine;     ///< Store each time used for refine timers
 
  bool isGeneric;     ///< define if the generic timer is used
  bool isNeighbours;  ///< define if the neighbour timer is used
  bool isPotential;   ///< define if the potential timer is used
  bool isGhost;       ///< define if the ghost timer is used
  bool isRefine;      ///< define if the refinement timer is used

  int number_of_step;  ///< number of step currently done

  // methods
  public :

  /// @brief constructor by default
  metricsDetails() : m_time_generic(0), m_time_neighbours(0), m_time_potential(0), m_time_ghost(0), m_time_refine(0), 
                     isGeneric(0), isNeighbours(0), isPotential(0), isGhost(0), isRefine(0), number_of_step(0)  {

    number_of_step = 0;

    // why not
    isGeneric    = false;
    isNeighbours = false;
    isPotential  = false;
    isGhost      = false;
    isRefine     = false;   

  }

  /// @brief accessor
  /// @return return booleen value
  bool isTimerNeighbours() {
    return isNeighbours;
  }

  /// @brief accessor
  /// @return return booleen value
  bool isTimerGeneric() {
    return isGeneric;
  }


  /// @brief accessor 
  /// @return return booleen value
  bool isTimerPotential() {
    return isPotential;
  }

  /// @brief accessor 
  /// @return return booleen value
  bool isTimerGhost() {
    return isGhost;
  }

  /// @brief accessor 
  /// @return return booleen value
  bool isTimerRefine() {
    return isRefine;
  }

  /// @brief Decide if the timer is chosen
  /// @param [in] b booleen
  /// @return void
  void useNeighbours(bool b) {
    isNeighbours=b;
    
  }

  /// @brief Decide if the timer is chosen
  /// @param [in] b booleen
  /// @return void
  void usePotential(bool b) {
    isPotential=b;
    
  }

  /// @brief Decide if the timer is chosen
  /// @param [in] b booleen
  /// @return void
  void useGhost(bool b) {
    isGhost=b;
    
  }

  /// @brief Decide if the timer is chosen
  /// @param [in] b booleen
  /// @return void
  void useRefine(bool b) {
    isRefine=b;
    
  }

  /// @brief Decide if the timer generic is chosen
  /// @param [in] b booleen
  /// @return void
  void useGeneric(bool b) {
    isGeneric=b;
    
  }

  /// @brief resize all arrays
  /// @param [in] n new size
  /// @return void
  void reserve(int n)
  {
      m_time_generic.reserve(n);
      m_time_neighbours.reserve(n);
      m_time_potential.reserve(n);
      m_time_ghost.reserve(n);
      m_time_refine.reserve(n);
  }
  /// @brief resize all arrays
  /// @param [in] n new size
  /// @return void
  void resize(int n)
  {
    if(n>number_of_step)
    {
      m_time_generic.resize(n,0);
      m_time_neighbours.resize(n,0);
      m_time_potential.resize(n,0);
      m_time_ghost.resize(n,0);
      m_time_refine.resize(n,0);
    }
    number_of_step = n;
  }

  /// @brief reset the data structur
  /// @return void
  void clear() {
    number_of_step=0;
  }

  /// @brief prepare the data structur for a new step
  /// @return void
  void resize_one() {
      resize(number_of_step+1);
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] time Type of the measured time
  /// @return void
  void fill_m_time_generic(double time) {
    m_time_generic[number_of_step-1]=time;
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] time Type of the measured time
  /// @return void
  void fill_m_time_neighbours(double time) {
    m_time_neighbours[number_of_step-1]=time;
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] time Type of the measured time
  /// @return void
  void fill_m_time_potential(double time) {
    m_time_potential[number_of_step-1]=time;
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] time Type of the measured time
  /// @return void
  void fill_m_time_ghost(double time) {
    m_time_ghost[number_of_step-1]=time;
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] time Type of the measured time
  /// @return void
  void fill_m_time_refine(double time) {
    m_time_refine[number_of_step-1]=time;
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] i iteration
  /// @return time of m_time_analytics for iteration i 
  double get_m_time_generic(int i) {
    return m_time_generic[i];
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] i iteration
  /// @return time of m_time_neighbours for iteration i 
  double get_m_time_neighbours(int i) {
    return m_time_neighbours[i];
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] i iteration
  /// @return time of m_time_potential for iteration i 
  double get_m_time_potential(int i) {
    return m_time_potential[i];
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] i iteration
  /// @return time of m_time_ghost for iteration i 
  double get_m_time_ghost(int i) {
    return m_time_ghost[i];
  }

  /// @brief Measure time on chosen chrono
  /// @param [in] i iteration
  /// @return time of m_time_refine for iteration i 
  double get_m_time_refine(int i) {
    return m_time_refine[i];
  }

  /// @brief Write timers in a data file  
  /// @param [in] name, name of the new data file
  /// @return void
  void writeFile(std::string name)
  {
    std::ofstream fichier(name, std::ios::out);
    fichier << "# iteration ";

    if(isGeneric)
      fichier << "| name " ;

    if(isNeighbours)
      fichier << "| neighbours " ;

    if(isPotential)
      fichier << "| potential " ;

    if(isGhost)
      fichier << "| ghost " ;

    if(isRefine)
      fichier << "| refine " ;

    fichier << std::endl;

    for(int i = 0 ; i < number_of_step ; i++)
    {
      fichier << i ;

      if(isGeneric)
      fichier << " " << get_m_time_generic(i);

      if(isNeighbours)
      fichier << " " <<  get_m_time_neighbours(i) ;

      if(isPotential)
      fichier << " " << get_m_time_potential(i);

      if(isGhost)
      fichier << " " << get_m_time_ghost(i) ;

      if(isRefine)
      fichier << " " << get_m_time_refine(i);

       fichier << std::endl;

    }

  }
 
  #include "mpi.h"

  // end of the simulation (i currently don't care of the performance, I just that be correct)
  void updateMPI(double &min, double &max, double &mean, double *tmp)
  {
    int count;
    MPI_Comm_size(MPI_COMM_WORLD, &count);

    double local_value = 0.;

    for(int i = 0; i<number_of_step ; i++)
      local_value += tmp[i];

    MPI_Reduce( &local_value, &min , 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce( &local_value, &max , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
    MPI_Reduce( &local_value, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 

    mean /= count;
  }


  /// @brief Print the time in day, hour, minute and second
  /// @param [in] a time value (double)
  /// @return void
  void print(double a)
  {

    int d = auxFloor<int>( int(a) / 86400.) ;
    a -= d*86400;
    int h = auxFloor<int>( int(a) / 3600.);
    a-=h* 3600;
    int m = auxFloor<int>( int(a) / 60.);
    a-=m*60;

    std::cout << std::setw(2) << d << " d. " << std::setw(2) << h << " h. " << std::setw(2) << m  << " m. "  << a << " s. ";
  }

  /// @brief Summary of the timers
  /// @param [in] total time of the simulation
  /// @return void
  void write_end(double total, bool print2)
  { 

    if(print2)
      std::cout << std::endl << std::endl << "#  MPI Optional Timers " << std::endl << std::endl;

    if(isNeighbours)
    {
      double min,max,mean;
      updateMPI(min,max,mean,m_time_neighbours.data());
      if(print2)
      {
        std::cout << "# Neighbours \n"; 
        std::cout << "#   min  : "; print(min); std::cout<< "\n";
        std::cout << "#   mean : "; print(mean); std::cout<< "\n";
        std::cout << "#   max  : "; print(max); std::cout<< "\n";
	std::cout << "# " << std::endl;
      }
    }

    if(isPotential)
    {
      double min,max,mean;
      updateMPI(min,max,mean,m_time_potential.data());
      if(print2) {
        std::cout << "# Potential \n";
        std::cout << "#   min  : "; print(min); std::cout<< "\n";
        std::cout << "#   mean : "; print(mean); std::cout<< "\n";
        std::cout << "#   max  : "; print(max); std::cout<< "\n";
	std::cout << "# " << std::endl;
      }
    }

    if(isGhost)
    {
      double min,max,mean;
      updateMPI(min,max,mean,m_time_ghost.data());
      if(print2) {
        std::cout << "# Ghost  \n";
        std::cout << "#  min  : "; print(min); std::cout<< "\n";
        std::cout << "#  mean : "; print(mean); std::cout<< "\n";
        std::cout << "#  max  : "; print(max); std::cout<< "\n";
	std::cout << "# " << std::endl;
      }
    }

    if(isRefine)
    {
      double min,max,mean;
      updateMPI(min,max,mean,m_time_refine.data());
      if(print2) {
        std::cout << "# Refine \n";
        std::cout << "#  min  : "; print(min); std::cout<< "\n";
        std::cout << "#  mean : "; print(mean); std::cout<< "\n";
        std::cout << "#  max  : "; print(max); std::cout<< "\n";
	std::cout << "# " << std::endl;
      }
    }
  }
};
