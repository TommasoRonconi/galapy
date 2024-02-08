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
	  
          a [shape=box, color="#61C0E8", style=filled, label="Install the library", href=""];
          b [shape=box, color="#61C0E8", style=filled, label="Download the database", href=""];
          c [shape=box, color="#B2C600", style=filled, label="Generate the parameter file", href=":ref:`param_file`", target="_blank" ];
          d [shape=box, color="#B2C600", style=filled, label="Run the fit", href=":ref:`photometric_fit`"];
          e [shape=box, color="#E78A61", style=filled, label="Results analysis", href=""];
	  
	  a1 [shape=box, color="#61C0E8", label="pip install galapy", fontname="Verbatim", fontsize=12];
	  b1 [shape=box, color="#61C0E8", label="galapy-download-database", fontname="Verbatim", fontsize=12];
	  c1 [shape=box, color="#B2C600", label="galapy-genparams", fontname="Verbatim", fontsize=12];
	  d1 [shape=box, color="#B2C600", label="galapy-fit", fontname="Verbatim", fontsize=12];

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
	  
	  { rank = same; "a"; "a1"; }
	  { rank = same; "b"; "b1"; }
	  { rank = same; "c"; "c1"; }
	  { rank = same; "d"; "d1"; }


          b -> c [ltail=cluster0,lhead=cluster1, minlen=2., style=dashed];

       }

..
    digraph G {
          size = "7.5,7.5";
          compound=true;
	  nodesep=.05;
	  fontname="Helvetica";
	  
          a [shape=box, color="#61C0E8", style=filled, label="Install the library", href=""];
          b [shape=box, color="#61C0E8", style=filled, label="Download the database", href=""];
          c [shape=box, color="#B2C600", style=filled, label="Generate the parameter file", href=":ref:`param_file`", target="_blank" ];
          d [shape=box, color="#B2C600", style=filled, label="Run the fit", href=":ref:`photometric_fit`"];
          e [shape=box, color="#E78A61", style=filled, label="Results analysis", href=""];

	  a0 [shape=point];
	  b0 [shape=point];
	  c0 [shape=point];
	  d0 [shape=point];
	  
	  a1 [shape=box, color="#61C0E8", style=dashed, label="pip install galapy", fontname="Verbatim"];
	  b1 [shape=box, color="#61C0E8", style=dashed, label="galapy-download-database", fontname="Verbatim"];
	  c1 [shape=box, color="#B2C600", style=dashed, label="galapy-genparams", fontname="Verbatim"];
	  d1 [shape=box, color="#B2C600", style=dashed, label="galapy-fit", fontname="Verbatim"];

	  subgraph cluster0 {
	     shape=plaintext;
             a -> b [maxlen=1.];
	     a -> a0 -> a1 [arrowhead=none];
	     b -> b0 -> b1 [arrowhead=none];
	     label="Preliminaries";
	     { rank = same; "a"; "a1"; "a0"; }
	     { rank = same; "b"; "b1"; "b0"; }
	  }
	  
          subgraph cluster1 {
             color=black;
	     c -> d -> e;
	     c -> c0 -> c1 [arrowhead=none, minsep=2];
	     d -> d0 -> d1 [arrowhead=none, minsep=2];
	     label="Workflow";
	     { rank = same; "c"; "c1"; "c0"; }
	     { rank = same; "d"; "d1"; "d0"; }
          }

          b -> c [ltail=cluster0,lhead=cluster1, minlen=2.];

       }
