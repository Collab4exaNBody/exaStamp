
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


