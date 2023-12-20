/// @file
/// @brief Reimplementation of some HerculeBase methods and some ParticleWriterHercule method that is here for unknown reasons


/// @brief [ECom] HerculeBase::initHerculeParticleDep
void HerculeBase::initHerculeParticleDep()
{
	if (!m_api.isNull()) return;
	/* initialisation de l'Api */
	m_api=HIc_Api("myApi","stdDif");

		HIc_Init_Standard_Services(m_api);
		HIc_Init_Standard_CEA_DAM_DIF(m_api, "");
		HIc_Init_Parallel_MPI_Services(m_api, comm->getCommunicator());

}


/// @brief [ECom] HerculeBase::openBase
void HerculeBase::openBase()
{
	if (!m_base.isNull()) return;

		m_base = HIc_Base(m_api, "myBase", "parallele");

	m_base.setItemConf("mode","write");
		m_base.setItemConf("write_dir",m_rep.c_str());
		m_base.setItemConf("read_dir",m_rep.c_str());
		m_base.setItemConf("bd_name","HDep");
		m_base.open();
}


/// @brief [ECom] HerculeBase::closeBase
void HerculeBase::closeBase()
{
	m_base.close();
	m_base=HIc_Base();
}


/// @brief [ECom] HerculeBase::openStep
void HerculeBase::openStep(uint _step)
{
	particles=buff;
		m_ctx = m_base.createCtxPar((double)_step, comm->getNumberOfNodes(), comm->getRank());

	m_ctx.setItemConf("mode","write");

	m_ctx.open();


}


/// @brief [ECom] ParticleWriterHercule::closeStep
void ParticleWriterHercule::closeStep()
{
	m_ctx.close();
    m_ctx=HIc_Ctx();
	particles=NULL;
}
