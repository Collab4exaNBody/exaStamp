/// @file
/// @brief Extension for input.hpp that defines and prints the input and output types

#ifndef _IO_H
#define _IO_H


/// @brief Enumeration of the output file types
enum OutputType {
	DAT, ///< DAT output file
	VTK, ///< VTK output file
	HERCULE_PARTICULE_DEP ///< Hercule output file
};


/// @brief Enumeration of the dump file types
enum DumpType {
	HERCULE_DUMP, ///< Hercule dump file
	LEGACY_DUMP, ///< Stamp-like dump file
	MOL_DUMP ///< Molecular file (same_format as molecular input)
};


/// @brief Enumeration of the initialization types
enum InitType {
	INIT_DEFAULT, ///< Initialization from a lattice
	INIT_STAMP_LEGACY_DUMP, ///< Initialization from a Stamp dump file
	INIT_HERCULE_DUMP, ///< Initialization from a Hercule dump file
	INIT_STAMPV4 ///< Initialization from a Stampv4
/*, INIT_HERCULE_DUMP_V1,INIT_HERCULE_DUMP_V2,INIT_HERCULE_DUMP_V3*/ };


/// @brief Stream insertion operator for an OutputType
/// @param [in,out] out Stream
/// @param [in] type OutputType to print
inline std::ostream& operator << (std::ostream& out, const OutputType& type){
	if (type == OutputType::DAT){
		out << "DAT";
	}else if (type == OutputType::VTK){
		out << "VTK";
	}else if (type == OutputType::HERCULE_PARTICULE_DEP){
		out << "hercule";
	}else{
		out << "???";
	}
	return out;
}


/// @brief Stream insertion operator for an DumpType
/// @param [in,out] out Stream
/// @param [in] type DumpType to print
inline std::ostream& operator << (std::ostream& out, const DumpType& type){
	if (type == DumpType::LEGACY_DUMP){
		out << "legacy";
	}else if (type == DumpType::HERCULE_DUMP){
		out << "hercule";
	}else{
		out << "???";
	}
	return out;
}


/// @brief Stream insertion operator for an InitType
/// @param [in,out] out Stream
/// @param [in] type InitType to print
inline std::ostream& operator << (std::ostream& out, const InitType& type){
	if (type == InitType::INIT_STAMP_LEGACY_DUMP){
		out << "legacy";
	}else if (type == InitType::INIT_HERCULE_DUMP){
		out << "hercule";
	}else if (type == InitType::INIT_STAMPV4){
		out << "stampv4";
	}else if (type == InitType::INIT_DEFAULT){
		out << "lattice";
	}else{
		out << "???";
	}
	return out;
}

#endif
