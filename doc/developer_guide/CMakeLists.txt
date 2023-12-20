# Build the developer's guide

set(LATEX_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
file(GLOB_RECURSE input_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.tex)
file(GLOB_RECURSE bib_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.bib)

# automatically convert odp files to pdf
file(GLOB_RECURSE odp_files ${CMAKE_CURRENT_SOURCE_DIR} *.odp)
set(GENERATED_IMAGES_DIR ${CMAKE_CURRENT_BINARY_DIR}/fig)
file(MAKE_DIRECTORY ${GENERATED_IMAGES_DIR})
unset(ODP_GENERATED_IMAGES)
foreach(odp_file ${odp_files})
	get_filename_component(basename ${odp_file} NAME)
	add_custom_command(OUTPUT ${GENERATED_IMAGES_DIR}/${basename} COMMAND soffice --headless --convert-to pdf ${odp_file} --outdir ${GENERATED_IMAGES_DIR})
	list(APPEND ODP_GENERATED_IMAGES ${GENERATED_IMAGES_DIR}/${basename})
endforeach()

ADD_LATEX_DOCUMENT(
	developer_guide.tex
	BIBFILES ${bib_files}
	INPUTS ${input_files}
	IMAGE_DIRS fig
	FORCE_PDF
	DEPENDS ${ODP_GENERATED_IMAGES}
	)

install(FILES ${LATEX_OUTPUT_PATH}/developer_guide.pdf DESTINATION doc)

