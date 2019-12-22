(* Paclet Info File *)

(* created 2019/11/11*)

Paclet[
  Name -> "TwoComponentSpinors",
  Version -> "0.1",
  MathematicaVersion -> "8+",
  Description -> "",
  Creator -> "",
  Extensions ->
  	{
  		{ (* Location of the Kernel *)
            "Kernel",
            Root -> FileNameJoin[{"TwoComponentSpinors", "Kernel"}],
            Context -> {"TwoComponentSpinors`"}
        },
        { (* Testing Files *)
            "Resource",
            Root->"Tests",
            Resources->{
                "WeylTrace.wl",
                "WeylLine.wl",
                "FermionSpinSum.wl"
            }
        },
  		{ (* Documentation files *)
            "Documentation",
            Language -> "English",
            MainPage -> "Guides/TwoComponentSpinors"
        }
  	}
]
