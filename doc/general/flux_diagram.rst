Quickstart
----------

.. graphviz::

    digraph process {
           a [shape=box, color="#B2C600", style=filled, label="first", href="http://google.com"];
           b [shape=box, color="#00B3C7", style=filled, label="second", href="https://www.wordreference.com/it/"];
           c [shape=box, color="#C71400", style=filled, label="third"];
           d [shape=box, color="#C700B3", style=filled, label="fourth"];
           a -> b -> c -> d;
        }
