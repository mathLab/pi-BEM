Parallel BEM Solver for deal.II
===============================

After you have compiled your application, you can run 

	make test

or
	
	ctest 

to start the testsuite.

Take a look at
https://www.dealii.org/developer/developers/testsuite.html for more
information on how to create tests and add categories of tests.

Notice to developers
====================

Before making a pull request, please make sure you run the script

    ./scripts/indent

from the base directory of this project, to ensure that no random 
white space changes are inserted in the repository.

The script requires Artistic Style Version 2.04 (astyle) to work 
properly.

Licence
=======

Please see the file ./LICENSE for details



