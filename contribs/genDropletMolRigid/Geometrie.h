/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/
// bulk
NumeroGeometrie = 1;
strcpy(ListeGeometries[NumeroGeometrie],"bulk");
NbParametresMailleGeometrie[NumeroGeometrie] = 0;
ListeTypeGeometrie[NumeroGeometrie] = GEO_BULK;

// aleatoire
NumeroGeometrie = 2;
strcpy(ListeGeometries[NumeroGeometrie],"aleatoire");
NbParametresMailleGeometrie[NumeroGeometrie] = 1;
ListeTypeGeometrie[NumeroGeometrie] = GEO_ALEATOIRE;

// Deux composes separés par un plan yz
NumeroGeometrie = 3;
strcpy(ListeGeometries[NumeroGeometrie],"bicomposant_selon_x");
NbParametresMailleGeometrie[NumeroGeometrie] = 10;
ListeTypeGeometrie[NumeroGeometrie] = GEO_BICOMPOSANT_SELON_X;

// Deux composes separés par un plan sinuosidal
NumeroGeometrie = 4;
strcpy(ListeGeometries[NumeroGeometrie],"rainure_sinus_un_materiau");
NbParametresMailleGeometrie[NumeroGeometrie] = 2;
ListeTypeGeometrie[NumeroGeometrie] = GEO_RAINURE_SINUS_UN_MATERIAU;

// Deux composes separés par un plan sinuosidal
NumeroGeometrie = 5;
strcpy(ListeGeometries[NumeroGeometrie],"bicomposant_selon_x_avec_rainure_sinus");
NbParametresMailleGeometrie[NumeroGeometrie] = 19;
ListeTypeGeometrie[NumeroGeometrie] = GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_SINUS;

// Deux composes separés par un plan triangulaire
NumeroGeometrie = 6;
strcpy(ListeGeometries[NumeroGeometrie],"bicomposant_selon_x_avec_rainure_triangulaire");
NbParametresMailleGeometrie[NumeroGeometrie] = 19;
ListeTypeGeometrie[NumeroGeometrie] = GEO_BICOMPOSANT_SELON_X_AVEC_RAINURE_TRIANGULAIRE;

// Sphere seule
NumeroGeometrie = 7;
strcpy(ListeGeometries[NumeroGeometrie],"sphere_seule");
NbParametresMailleGeometrie[NumeroGeometrie] = 4;
ListeTypeGeometrie[NumeroGeometrie] = GEO_SPHERE_SEULE;


