generate_default_species:
  - message: "Initialize H2O molecule"
  - species:
      verbose: true
      species:
        - H: { mass: 1.007 Da , z: 1 , charge: 0.5564 e-  , molecule: H2O }
        - O: { mass: 16 Da    , z: 8 , charge: 0.0 e-     , molecule: H2O }
        - M: { mass: 0.000 Da , z: 0 , charge: -1.1128 e- , molecule: H2O }
        - TIP4P:
            charge: 0.0 e-
            z: 0
            mass: 0.0 Da
            rigid_molecule:
              - O: [  0.00000000 ang ,  0.00000000 ang ,  0.00000000 ang ]
              - H: [  0.75695033 ang ,  0.58588228 ang ,  0.00000000 ang ]
              - H: [ -0.75695033 ang ,  0.58588228 ang ,  0.00000000 ang ]
              - M: [  0.00000000 ang ,  0.15460000 ang ,  0.00000000 ang ]

