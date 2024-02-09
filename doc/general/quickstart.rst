Quickstart
----------

.. graphviz::
   :align: center
   :caption: Flux diagram summarising the main steps required to run GalaPy.
	     Filled boxes can be clicked to be re-directed to the relative documentation
	     page with in-detail explanation of the step.

    digraph G {
          size = "7.5,7.5";
          compound=true;
	  nodesep=.05;
	  fontname="Helvetica";
	  newrank=true;
	  
          a [shape=box, color="#61C0E8", style=filled, label="Install the library", href="https://galapy.readthedocs.io/en/latest/general/install_guide.html", target="_blank"];
          b [shape=box, color="#61C0E8", style=filled, label="Download the database", href="https://galapy.readthedocs.io/en/latest/general/install_guide.html#after-install", target="_blank"];
          c [shape=box, color="#B2C600", style=filled, label="Generate the parameter file", href="https://galapy.readthedocs.io/en/latest/guides/parameter_file.html", target="_blank" ];
          d [shape=box, color="#B2C600", style=filled, label="Run the fit", href="https://galapy.readthedocs.io/en/latest/guides/photometric_fit.html", target="_blank"];
          e [shape=box, color="#E78A61", style=filled, label="Results analysis", href="https://galapy.readthedocs.io/en/latest/notebooks/results_analysis.html", target="_blank"];
	  
	  install_command [shape=box, color="#61C0E8", label="pip install galapy", fontname="Verbatim", fontsize=12];
	  download_command [shape=box, color="#61C0E8", label="galapy-download-database", fontname="Verbatim", fontsize=12];
	  genparams_command [shape=box, color="#B2C600", label="galapy-genparams", fontname="Verbatim", fontsize=12];
	  fit_command [shape=box, color="#B2C600", label="galapy-fit", fontname="Verbatim", fontsize=12];

	  subgraph cluster0 {
	     color="black";
	     style=dashed
             a -> b [maxlen=1.];
	     label="GalaPy Library Set-Up";
	  }
	  
          subgraph cluster1 {
             color=black;
	     c -> d -> e;
	     label="GalaPy Workflow";
          }
	  
	  { rank = same; "a"; "install_command"; }
	  { rank = same; "b"; "download_command"; }
	  { rank = same; "c"; "genparams_command"; }
	  { rank = same; "d"; "fit_command"; }


          b -> c [ltail=cluster0,lhead=cluster1, minlen=2., style=dashed];

       }
