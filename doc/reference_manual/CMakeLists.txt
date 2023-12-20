# Build reference manual

set(LATEX_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
file(GLOB_RECURSE input_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.tex)
file(GLOB_RECURSE bib_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.bib)

ADD_LATEX_DOCUMENT(
	reference_manual.tex
	BIBFILES ${bib_files}
	INPUTS ${input_files}
	IMAGE_DIRS fig
	FORCE_PDF
	)

install(FILES ${LATEX_OUTPUT_PATH}/reference_manual.pdf DESTINATION doc)

