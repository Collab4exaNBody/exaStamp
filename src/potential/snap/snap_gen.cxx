#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cstdio>
#include <string>
#include <algorithm>
#include <cmath>
#include <cassert>

std::string output_dir = ".";

struct UExpr
{
  int idx = -1;
  std::string expr;
  std::set<int> deps;
  
  inline int conjsrc() const
  {
    if( deps.size() != 1 ) return -1;
    int src = * deps.begin();
    char buf[256];
    snprintf(buf,256,"conj(U_%03d)",src);
    if(expr==buf) return src;
    else return -1;
  }
  
  inline int negconjsrc() const
  {
    if( deps.size() != 1 ) return -1;
    int src = * deps.begin();
    char buf[256];
    snprintf(buf,256,"-conj(U_%03d)",src);
    if(expr==buf) return src;
    else return -1;
  }
};

void sreplace(std::string& s , const std::string& from , const std::string& to)
{
  auto p = s.find(from);
  while( p != std::string::npos )
  {
    s.replace( p , from.length() , to );
    p = s.find(from,p+to.length());
  }
}

struct SNA_ZINDICES {
  int16_t j1;
  int16_t j2;
  int16_t j;
  int16_t ma1min;
  int16_t ma2max;
  int16_t mb1min;
  int16_t mb2max;
  int16_t na;
  int16_t nb;
  int16_t jju;
};


void uarray_code_gen(int jmax , bool verbose = false )
{
  const int twojmax = 2 * jmax;

  std::vector<int> idxu_block( twojmax+1 , -1 );
  int idxu_count = 0;
  for (int j = 0; j <= twojmax; j++) {
    idxu_block[j] = idxu_count;
    for (int mb = 0; mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  const int idxu_max = idxu_count;

  char buf[256];
  std::vector< std::vector<UExpr> > ulist( twojmax+1 );
  std::vector<UExpr> prog;

  char rpqvar[256];
  std::map<std::string,double> rpq_map;
  auto rootpq = [&rpqvar,&rpq_map](int p, int q) -> const char*
  {
    int maxdiv=1;
    for(int d=2;d<=std::min(p,q);d++)
    {
      if( d*(p/d) == p && d*(q/d) == q ) maxdiv = d;
    }
    p /= maxdiv;
    q /= maxdiv;

    const double rpq = std::sqrt(static_cast<double>(p)/q);
    if( rpq==1.0 ) return "1.0";
    else if ( rpq==0.5 ) return "0.5";
    else if ( rpq==2.0 ) return "2.0";
    else
    {
      snprintf(rpqvar,256,"rpq_%d_%d",p,q);
      rpq_map[rpqvar] = rpq;
      return rpqvar;
    }
  };

  // special post processing case for derivatives, fot those left elements :
  std::set<int> dU_postProcess;
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        dU_postProcess.insert(jju);
        jju++;
      }
  }

  // dU and Y terms used for dE/dr (force) calculation
  std::set<int> dU_Y_unscaled;
  std::set<int> dU_Y_scaled;
  std::set<int> dU_Y_unused;
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int mb = 0; 2*mb < j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        dU_Y_unscaled.insert(jju);
        jju++;
      } //end loop over ma mb

    // For j even, handle middle column
    if (j%2 == 0) {
      int mb = j/2;
      for (int ma = 0; ma < mb; ma++) {
        dU_Y_unscaled.insert(jju);
        jju++;
      }
      dU_Y_scaled.insert(jju);
      // jju++;
    } // end if jeven
  } // end loop over j

  int force_used_jju = 0;
  for(int jju=0;jju<idxu_max;jju++)
  {
    if( dU_Y_unscaled.find(jju)!=dU_Y_unscaled.end() || dU_Y_scaled.find(jju)!=dU_Y_scaled.end() )
    {
      force_used_jju ++;
    }
    else
    {
      dU_Y_unused.insert( jju );
    }
  }
  
  if( verbose )
  {
    std::cout<<"force uses "<<force_used_jju <<" / "<<idxu_max <<" dU*Y terms"<<std::endl;
    std::cout<<"2x scaled terms :";
    for(const auto jju:dU_Y_unscaled) std::cout<<" "<<jju;
    std::cout<<std::endl;
    std::cout<<"1x scaled terms :";
    for(const auto jju:dU_Y_scaled) std::cout<<" "<<jju;
    std::cout<<std::endl;
    std::cout<<"unused terms :";
    for(const auto jju:dU_Y_unused) std::cout<<" "<<jju;
    std::cout<<std::endl;
  }


  // Y terms dependencies on UTot terms
  int idxz_max = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_max++;

  std::vector<SNA_ZINDICES> IDXZ(idxz_max);

  int idxz_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            IDXZ[idxz_count].j1 = j1;
            IDXZ[idxz_count].j2 = j2;
            IDXZ[idxz_count].j = j;
            IDXZ[idxz_count].ma1min = std::max(0, (2 * ma - j - j2 + j1) / 2);
            IDXZ[idxz_count].ma2max = (2 * ma - j - (2 * IDXZ[idxz_count].ma1min - j1) + j2) / 2;
            IDXZ[idxz_count].na = std::min(j1, (2 * ma - j + j2 + j1) / 2) - IDXZ[idxz_count].ma1min + 1;
            IDXZ[idxz_count].mb1min = std::max(0, (2 * mb - j - j2 + j1) / 2);
            IDXZ[idxz_count].mb2max = (2 * mb - j - (2 * IDXZ[idxz_count].mb1min - j1) + j2) / 2;
            IDXZ[idxz_count].nb = std::min(j1, (2 * mb - j + j2 + j1) / 2) - IDXZ[idxz_count].mb1min + 1;
            const int jju = idxu_block[j] + (j+1)*mb + ma;
            IDXZ[idxz_count].jju = jju;
            idxz_count++;
          }
      }

  std::vector< std::set<int> > Y_Utot_deps( idxu_max );
  std::set<int> Y_computed;
  for (int jjz = 0; jjz < idxz_max; jjz++)
  {
    const int j1 = IDXZ[jjz].j1;
    const int j2 = IDXZ[jjz].j2;
    const int j = IDXZ[jjz].j;
    const int ma1min = IDXZ[jjz].ma1min;
    const int ma2max = IDXZ[jjz].ma2max;
    const int na = IDXZ[jjz].na;
    const int mb1min = IDXZ[jjz].mb1min;
    const int mb2max = IDXZ[jjz].mb2max;
    const int nb = IDXZ[jjz].nb;
    const int jju = IDXZ[jjz].jju;
    Y_computed.insert( jju );
    int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
    int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
    for (int ib = 0; ib < nb; ib++)
    {
      int ma1 = ma1min;
      int ma2 = ma2max;
      for (int ia = 0; ia < na; ia++)
      {
        Y_Utot_deps[jju].insert( jju1+ma1 );
        Y_Utot_deps[jju].insert( jju2+ma2  ); 
        ma1++;
        ma2--;
      } // end loop over ia
      jju1 += j1 + 1;
      jju2 -= j2 + 1;
    } // end loop over ib
  } // end loop over jjz
  
  int Y_jju_count = 0;
  std::map<int,int> Y_jju_map;
  std::set<int> Y_wasted;
  for(auto idx:Y_computed)
  {
    if( dU_Y_unused.find(idx) != dU_Y_unused.end() ) Y_wasted.insert(idx);
    else Y_jju_map[ idx ] = Y_jju_count ++;
  }

  auto y_jju_mapped_index = [&Y_jju_map](int jju) -> int
  {
    if( Y_jju_map.find(jju) == Y_jju_map.end() ) return -1;
    else return Y_jju_map[jju];
  };

  std::vector<bool> UTot_used_for_Y( idxu_max , false );
  for(int jju=0;jju<idxu_max;jju++)
  {
    UTot_used_for_Y[jju] = UTot_used_for_Y[jju] || ( Y_Utot_deps[jju].find(jju) != Y_Utot_deps[jju].end() );
  }
  
  {
    std::string ymapfilename = output_dir + "/ymap_jmax"+std::to_string(jmax)+".hxx";
    if( std::ifstream( ymapfilename ).good() )
    {
      std::cout<<"skip "<<ymapfilename<<std::endl;
    }
    else
    {
      std::cout<<"generate "<<ymapfilename<<" ..."<<std::endl;
      std::ofstream ymapfile( ymapfilename );
      ymapfile << "static constexpr int Y_jmax"<<jmax<<"_jju_count = "<<Y_jju_count<<";"<<std::endl;
      ymapfile << "static constexpr int Y_jmax"<<jmax<<"_jju_map["<<idxu_max+1<<"] = {"<<std::endl;
      for(int jju=0;jju<idxu_max;jju++) ymapfile << "\t"<< y_jju_mapped_index(jju) <<" ," << std::endl;
      ymapfile << "\t-1 };" << std::endl;
    }
  }
    
  if( verbose )
  {
    std::cout << "computed Y terms :";
    for(auto idx:Y_computed) std::cout<<" "<<idx;
    std::cout<<std::endl;

    std::cout << "UTot terms used by Y computation :";
    int utot_count=0;
    for(int jju=0;jju<idxu_max;jju++) if( UTot_used_for_Y[jju] ) { std::cout<<" "<<jju; ++utot_count; }
    std::cout<<std::endl<<"Y uses "<<utot_count<<" / "<<idxu_max<<" terms of UTot"<<std::endl;
    /*
    for(int jju=0;jju<idxu_max;jju++)
    {
      std::cout<<"Y("<<jju<<") depends on UTot indices :";
      for(auto idx:Y_Utot_deps[jju]) std::cout<<" "<<idx;
      std::cout<<std::endl;
    }
    */
    std::cout << "Y terms computed but useless for force :";
    for(auto idx:Y_wasted) std::cout<<" "<<idx;
    std::cout<<std::endl;
  }


  /*****************************
   * Expression list generation
   *****************************/

//  printf("U_%03d = (1,0)\n",0);
  
  ulist[0].push_back( { 0 , "U_UNIT" , {} } );

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];
    //printf("// j=%d : idxu_block = [%d;%d[\n",j,jjup,jju);
    for (int mb = 0; 2*mb <= j; mb++) {
      //printf("U_%03d = (0,0)\n",jju);
      ulist[j].push_back( { jju, "U_ZERO" , {} } );
      for (int ma = 0; ma < j; ma++) {
        // U(jju) += rootpq(j-ma,j-mb) * conj(a) * U(jjup)
//        snprintf(buf,256,"U_%03d + %s * ( conj(a) * U_%03d )",jju, rootpq(j-ma,j-mb) ,jjup );
        snprintf(buf,256,"U_%03d + %s * U_BLEND(a,U_%03d)",jju, rootpq(j-ma,j-mb) ,jjup );
        ulist[j].push_back( { jju , buf , {jju,jjup} } );
        // U(jju+1) = -rootpq(ma+1,j-mb) * conj(b) * U(jjup)
        snprintf(buf,256,"-%s * U_BLEND(b,U_%03d)", rootpq(ma+1,j-mb) ,jjup );
        ulist[j].push_back( { jju+1 , buf , {jjup} } );
        jju++;
        jjup++;
      }
      jju++;
    }

    jju = idxu_block[j];
    jjup = jju+(j+1)*(j+1)-1;
    int mbpar = 1;
    for (int mb = 0; 2*mb <= j; mb++) {
      int mapar = mbpar;
      for (int ma = 0; ma <= j; ma++) {
//      U(jjup) = +/- conj( U(jju) )
        snprintf(buf,256,"%sconj(U_%03d)",((mapar==1)?"":"-"),jju);
        ulist[j].push_back( { jjup , buf , {jju} } );
        mapar = -mapar;
        jju++;
        jjup--;
      }
      mbpar = -mbpar;
    }
  }
  
  int block=0;
  for(auto & ul : ulist)
  {
    int bs = idxu_block[block];
    int be = block<twojmax ? idxu_block[block+1] : idxu_max ;
    if(verbose) std::cout<<"--- J="<<block<<" block ["<<bs<<";"<<be<<"[ ---"<<std::endl;
    block++;

    if(verbose) for(const auto & e : ul )
    {
      snprintf(buf,256,"U_%03d",e.idx);
      std::cout << buf << " = "<<e.expr<<"; // depends on";
      for(const auto& d : e.deps) std::cout<<" "<<d;
      std::cout<<std::endl;
    }

    if(verbose) std::cout<<"--- substitution ---"<<std::endl;

    const int n = ul.size();
    for(int i=0;i<n;i++)
    {
      bool supstitute = false;
      if( ul[i].deps.empty() ) supstitute = true;
      else
      {
        for(int j=i+1;j<n;j++)
        {
          if( ul[j].idx == ul[i].idx )
          {
            if( ul[j].deps.find(ul[i].idx) != ul[j].deps.end() ) supstitute = true;
            break;
          }
        }
      }
      if( supstitute )
      {
        snprintf(buf,256,"U_%03d",ul[i].idx);
        std::string Ustr = buf;
        int replace_up_to = ul.size() - 1;
        for(int j=i+1;j<n;j++) if( ul[j].idx==ul[i].idx ) { replace_up_to=j; break; }
//        std::cout << "replace U_"<<ul[i].idx<<" with "<<ul[i].expr<<" in range ["<<i+1<<";"<<replace_up_to<<"]"<<std::endl;
        for(int j=i+1;j<=replace_up_to;j++)
        {
          if( ul[j].deps.find(ul[i].idx) != ul[j].deps.end() )
          {
            auto p = ul[j].expr.find(Ustr);
            std::string oexpr = ul[j].expr;
            ul[j].expr.replace(p,Ustr.length(),ul[i].expr);
            ul[j].deps.erase( ul[i].idx );
            ul[j].deps.insert( ul[i].deps.begin() , ul[i].deps.end() );
            if(verbose) std::cout<<"expr #"<<j<<" for U_"<<ul[j].idx<<" : "<<oexpr<<" -> "<<ul[j].expr<<std::endl;
          }
        }
      }
    }

    if(verbose) for(const auto & e : ul )
    {
      snprintf(buf,256,"U_%03d",e.idx);
      std::cout << buf << " = "<<e.expr<<"; // depends on";
      for(const auto& d : e.deps) std::cout<<" "<<d;
      std::cout<<std::endl;
    }

    if(verbose) std::cout<<"--- deletion ---"<<std::endl;
    for(int i=0;i<n;i++)
    {
      bool deletable = false;
      for(int j=i+1;j<n;j++)
      {
        if( ul[j].deps.find(ul[i].idx) != ul[j].deps.end() ) break;
        if( ul[j].idx == ul[i].idx ) { deletable=true; break; }
      }
      if( deletable ) ul[i].idx = - ul[i].idx;
    }
    if(verbose) for(const auto & e : ul )
    {
      if(e.idx>=0) 
      {
        snprintf(buf,256,"U_%03d",e.idx);
        std::cout << buf << " = "<<e.expr<<";";
        if(!e.deps.empty()) std::cout<<" // depends on";
        for(const auto& d : e.deps) std::cout<<" "<<d;
        std::cout<<std::endl;
      }
    }

    if(verbose) std::cout<<"--- simplification ---"<<std::endl;
    for(int i=1;i<n;i++)
    {
      if( ul[i].idx >= 0 )
      {
        int src = ul[i].conjsrc();
        if( src != -1 )
        {
          int srcdef = -1;
          for(int j=i-1;j>=0;j--) if( ul[j].idx == src ) { srcdef=j; break; }
          if( srcdef != -1 )
          {
            if( ul[srcdef].conjsrc() == ul[i].idx )
            {
              if(verbose) std::cout<<"U_"<<ul[i].idx<<" is conjugate of U_"<<src<<" which is conjugate of U_"<<ul[i].idx<<std::endl;
              ul[i].idx = - ul[i].idx;
            }
          }          
        }
        else
        {
          src = ul[i].negconjsrc();
          if( src != -1 )
          {
            int srcdef = -1;
            for(int j=i-1;j>=0;j--) if( ul[j].idx == src ) { srcdef=j; break; }
            if( srcdef != -1 )
            {
              if( ul[srcdef].negconjsrc() == ul[i].idx )
              {
                if(verbose) std::cout<<"U_"<<ul[i].idx<<" is negative conjugate of U_"<<src<<" which is negative conjugate of U_"<<ul[i].idx<<std::endl;          
                ul[i].idx = - ul[i].idx;
              }
            }          
          }
        }
      }
    }

    // check if reorder is possible
    bool can_reorder = true;
    for(const auto & e : ul )
    {
      if(e.idx>=0) 
      {
        for(auto d:e.deps) if(d>=e.idx) can_reorder=false;
      }
    }
    if(verbose) std::cout<<"reorder="<<std::boolalpha<<can_reorder<<std::endl;
    if( can_reorder ) std::stable_sort( ul.begin() , ul.end() , [](const UExpr& a, const UExpr& b)->bool { return a.idx < b.idx; } );
    
    for(int i=0;i<n;i++)
    {
      if(ul[i].idx>=0) prog.push_back( ul[i] );
    }    
  }

  auto & ul = prog;
  int n = ul.size();  

  // final pass : compute variable life cycle end points and output code
  std::vector<int> last_use_pc( idxu_max , -1 );
  std::vector<int> first_set_pc( idxu_max , -1 );
  std::vector<int> last_set_pc( idxu_max , -1 );
  for(int i=0;i<n;i++)
  {
    if(ul[i].idx>=0)
    {
      last_use_pc.at(ul[i].idx)=i;
      last_set_pc.at(ul[i].idx)=i;
      for(auto d:ul[i].deps) last_use_pc.at(d)=i;
      if( first_set_pc.at(ul[i].idx) == -1 ) first_set_pc.at(ul[i].idx) = i;
    }
  }
  
  std::vector< std::set<int> > last_use_at( n );
  std::vector< std::set<int> > first_set_at( n );
  std::vector< std::set<int> > last_set_at( n );
  for(int jju=0;jju<idxu_max;jju++)
  {
    int lastuse = last_use_pc.at(jju);
    int firstset = first_set_pc.at(jju);
    int lastset = last_set_pc.at(jju);

    if( lastuse != -1 ) last_use_at.at(lastuse).insert( jju );
    if( firstset != -1 ) first_set_at.at(firstset).insert( jju );
    if( lastset != -1 ) last_set_at.at(lastset).insert( jju );
  }



  /**************
   * Code output 
   **************/
  std::string filename = output_dir + "/compute_ui_jmax" + std::to_string(jmax) + ".hxx";
  if( std::ifstream( filename ).good() )
  {
    std::cout<<"skip "<<filename<<std::endl;
    return ;
  }
  std::cout << "generate "<<filename<<" ..."<<std::endl;
  std::ofstream fout( filename );
  
  std::vector<std::string> macros = { "CONST_DECLARE(var,value)"
                                    , "U_BLEND(c,var)"
                                    , "U_BLEND_R(c,var)"
                                    , "U_BLEND_I(c,var)"
                                    , "U_FINAL(var,i)"
                                    , "U_STORE(var,i)"
                                    , "U_ASSIGN(var,expr)"
                                    , "U_ASSIGN_R(var,expr)"
                                    , "U_ASSIGN_I(var,expr)"
                                    , "U_DECLARE(var)"
                                    , "dU_DECLARE(var)"
                                    , "dU_ASSIGN(var,expr)"
                                    , "dU_ASSIGN_X_R(var,expr)"
                                    , "dU_ASSIGN_X_I(var,expr)"
                                    , "dU_ASSIGN_Y_R(var,expr)"
                                    , "dU_ASSIGN_Y_I(var,expr)"
                                    , "dU_ASSIGN_Z_R(var,expr)"
                                    , "dU_ASSIGN_Z_I(var,expr)"
                                    , "dU_BLEND(c,var)"
                                    , "dU_BLEND_X_R(c,var)"
                                    , "dU_BLEND_X_I(c,var)"
                                    , "dU_BLEND_Y_R(c,var)"
                                    , "dU_BLEND_Y_I(c,var)"
                                    , "dU_BLEND_Z_R(c,var)"
                                    , "dU_BLEND_Z_I(c,var)"
                                    , "dU_FINAL(var,i)"
                                    , "dU_STORE(var,i,j)"
                                    , "dU_POSTPROCESS(var,i)"
                                    , "BAKE_U_BLEND(c,var)"
                                    };

  std::map< std::string , std::string > macro_variants = {
      { "U_STORE_FSKIP"        , "U_STORE"       }
    , { "dU_STORE_FHALF"       , "dU_STORE"      } 
    , { "dU_STORE_FSKIP"       , "dU_STORE"      } 
    , { "dU_POSTPROCESS_FSKIP" , "dU_POSTPROCESS"} };

  for(const auto& m : macros)
  {
    std::string name = m.substr( 0, m.find('(') );
    fout<<"#ifndef "<<name<<std::endl;
    fout<<"#define "<<m<<" /**/"<<std::endl;
    fout<<"#endif"<<std::endl;
  }
  for(const auto& m : macro_variants)
  {
    fout<<"#ifndef "<<m.first<<std::endl;
    fout<<"#define "<<m.first<<" "<<m.second<<std::endl;
    fout<<"#endif"<<std::endl;
  }

  for(const auto& rpq : rpq_map)
  {
    fout<<"CONST_DECLARE("<<rpq.first<<","<<std::hexfloat<<rpq.second<<");"<<std::endl;
  }

  std::set<std::string> precomputed_blend;
  auto precompute_u_blend = [&](const std::string& expr)
  {
    static const std::string BSTR = "U_BLEND(";
    auto s = expr.find(BSTR);
    while( s != std::string::npos )
    {
      auto e = expr.find(")",s+BSTR.length());
      assert( e != std::string::npos );
      std::string b = expr.substr(s,e-s+1);
      if( precomputed_blend.find(b) == precomputed_blend.end() )
      {
        fout<<"BAKE_"<<b<<";"<<std::endl;
        precomputed_blend.insert(b);
      }
      s = expr.find(BSTR,e+1);
    }
  };

  bool multiple_set = false;
  for(int i=0;i<n;i++) if( ul[i].idx >= 0 )
  {
    snprintf(buf,256,"U_%03d",ul[i].idx);
    if( first_set_at[i].find(ul[i].idx) != first_set_at[i].end() )
    {
      fout<<"U_DECLARE("<<buf<<");"<<std::endl;
      fout<<"dU_DECLARE("<<buf<<");"<<std::endl;
    }
    else fout<<"// overwrite of "<<buf<<std::endl;
    std::string u_expr = ul[i].expr;
    sreplace(u_expr,"U_ZERO + ","");
    sreplace(u_expr,"1.0 * ","");
    sreplace(u_expr," * ","*");
    sreplace(u_expr," + ","+");
    precompute_u_blend(u_expr);

    fout << "U_ASSIGN("<<buf<<","<<u_expr<<");"<<std::endl;
    
    std::string u_r_expr = u_expr;
    sreplace(u_r_expr,"U_ZERO","U_ZERO_R");
    sreplace(u_r_expr,"U_UNIT","U_UNIT_R");
    sreplace(u_r_expr,"conj(U_BLEND","CONJ_R(U_BLEND"); 
    sreplace(u_r_expr,"conj(U_","CONJ_U_R(U_"); 
    sreplace(u_r_expr,"conj(","CONJ_R("); 
    sreplace(u_r_expr,"U_BLEND(","U_BLEND_R("); 
    fout << "U_ASSIGN_R("<<buf<<","<<u_r_expr<<");"<<std::endl;

    std::string u_i_expr = u_expr;
    sreplace(u_i_expr,"U_ZERO","U_ZERO_I");
    sreplace(u_i_expr,"U_UNIT","U_UNIT_I");
    sreplace(u_i_expr,"conj(U_BLEND","CONJ_I(U_BLEND"); 
    sreplace(u_i_expr,"conj(U_","CONJ_U_I(U_"); 
    sreplace(u_i_expr,"conj(","CONJ_I("); 
    sreplace(u_i_expr,"U_BLEND(","U_BLEND_I("); 
//    sreplace(u_i_expr,"-CONJ_U_I(","U_I_ID("); 
//    sreplace(u_i_expr,"-CONJ_I(","U_ID("); 
    fout << "U_ASSIGN_I("<<buf<<","<<u_i_expr<<");"<<std::endl;

    std::string du_expr = u_expr;
    sreplace(du_expr,"U_ZERO","dU_ZERO");
    sreplace(du_expr,"U_UNIT","dU_ZERO");
    sreplace(du_expr,"U_BLEND","dU_BLEND"); 
    sreplace(du_expr,"conj(U_","CONJC3D(dU_"); 
    sreplace(du_expr,"conj(","CONJC3D("); 
    fout << "dU_ASSIGN("<<buf<<","<<du_expr<<");"<<std::endl;
    
    std::string du_x_expr = du_expr;
    sreplace(du_x_expr,"dU_ZERO","dU_ZERO_X_R");
    sreplace(du_x_expr,"dU_UNIT","dU_UNIT_X_R");
    sreplace(du_x_expr,"CONJC3D(dU_BLEND","CONJC3D_X_R(dU_BLEND"); 
    sreplace(du_x_expr,"CONJC3D(dU_","CONJC3D_dU_X_R(dU_"); 
    sreplace(du_x_expr,"CONJC3D(","CONJC3D_X_R("); 
    sreplace(du_x_expr,"dU_BLEND(","dU_BLEND_X_R("); 
    fout << "dU_ASSIGN_X_R("<<buf<<","<<du_x_expr<<");"<<std::endl;
    sreplace(du_x_expr,"_X_R","_X_I");
    fout << "dU_ASSIGN_X_I("<<buf<<","<<du_x_expr<<");"<<std::endl;
    sreplace(du_x_expr,"_X_I","_Y_R");
    fout << "dU_ASSIGN_Y_R("<<buf<<","<<du_x_expr<<");"<<std::endl;
    sreplace(du_x_expr,"_Y_R","_Y_I");
    fout << "dU_ASSIGN_Y_I("<<buf<<","<<du_x_expr<<");"<<std::endl;
    sreplace(du_x_expr,"_Y_I","_Z_R");
    fout << "dU_ASSIGN_Z_R("<<buf<<","<<du_x_expr<<");"<<std::endl;
    sreplace(du_x_expr,"_Z_R","_Z_I");
    fout << "dU_ASSIGN_Z_I("<<buf<<","<<du_x_expr<<");"<<std::endl;
    
    if( last_set_at[i].find(ul[i].idx) != last_set_at[i].end() )
    {
      const char * s = ( dU_Y_unused.find(ul[i].idx) != dU_Y_unused.end() ) ? "_FSKIP" : "";
      fout<<"U_STORE"<<s<<"("<<buf<<","<<ul[i].idx<<");"<<std::endl;
      if( dU_postProcess.find(ul[i].idx) == dU_postProcess.end() )
      {
        const char * sdu = ( dU_Y_scaled.find(ul[i].idx) != dU_Y_scaled.end() ) ? "_FHALF" : s;
        fout<<"dU_STORE"<<sdu<<"("<<buf<<","<<ul[i].idx<<","<<y_jju_mapped_index(ul[i].idx)<<");"<<std::endl;
      }
    }
    for(int dvar:last_use_at[i])
    {
      snprintf(buf,256,"U_%03d",dvar);
      if( dU_postProcess.find(dvar) != dU_postProcess.end() )
      {
        const char * s = ( dU_Y_unused.find(dvar) != dU_Y_unused.end() ) ? "_FSKIP" : "";      
        const char * sdu = ( dU_Y_scaled.find(dvar) != dU_Y_scaled.end() ) ? "_FHALF" : s;
        fout<<"dU_POSTPROCESS"<<s<<"("<<buf<<","<<dvar<<");"<<std::endl;
        fout<<"dU_STORE"<<sdu<<"("<<buf<<","<<dvar<<","<<y_jju_mapped_index(dvar)<<");"<<std::endl;
      }
      fout<<"U_FINAL("<<buf<<","<<dvar<<");"<<std::endl;
      fout<<"dU_FINAL("<<buf<<","<<dvar<<");"<<std::endl;
    }
  }

  fout<<"#ifndef SNAP_AUTOGEN_NO_UNDEF"<<std::endl;
  for(const auto& m : macros)
  {
    std::string name = m.substr( 0, m.find('(') );
    fout<<"#undef "<<name<<std::endl;
  }
  for(const auto& m : macro_variants)
  {
    fout<<"#undef "<<m.first<<std::endl;
  }
  fout<<"#endif"<<std::endl;
  
}


int main(int argc, char* argv[] )
{
  bool verbose = false;
  if( argc >= 2 ) output_dir = argv[1];
  if( argc == 3 ) verbose = ( std::string(argv[2]) == "-v" );
  for(int jmax=2;jmax<=4;jmax++) uarray_code_gen(jmax,verbose);
  return 0;
}

