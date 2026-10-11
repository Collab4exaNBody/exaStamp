/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/


#pragma once

#include <string>
#include <memory>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>

namespace exaStamp
{
  // Old operator name kept for existing input files : prints a warning (once), then builds the operator registered as
  // new_name with the same YAML configuration.
  inline void register_deprecated_operator_alias( const std::string& old_name, const std::string& new_name )
  {
    using namespace onika::scg;
    auto warned = std::make_shared<bool>( false );
    OperatorNodeFactory::instance()->register_factory( old_name ,
      [old_name,new_name,warned]( const YAML::Node& node, const OperatorNodeFlavor& flavor ) -> std::shared_ptr<OperatorNode>
      {
        if( ! *warned ) { onika::lout << "Warning: operator '"<< old_name <<"' is deprecated, use '"<< new_name <<"' instead" << std::endl; *warned = true; }
        return OperatorNodeFactory::instance()->make_operator( new_name , node , flavor );
      } );
  }

  // Operator name that no longer exists, and whose replacement does not take the same parameters :
  // building it is a fatal error explaining what to use instead.
  inline void register_removed_operator( const std::string& old_name, const std::string& hint )
  {
    using namespace onika::scg;
    OperatorNodeFactory::instance()->register_factory( old_name ,
      [old_name,hint]( const YAML::Node&, const OperatorNodeFlavor& ) -> std::shared_ptr<OperatorNode>
      {
        onika::fatal_error() << "operator '"<< old_name <<"' has been removed. "<< hint << std::endl;
        return nullptr;
      } );
  }
}
