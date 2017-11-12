import os

from Quadrature_Analysis import main

def test_main_function_with_correct_arguments():
    assert main('f1', '1', '1', '10') == None

def test_errors_produced_by_f1():
    errors_file_path = os.path.join('errors', 'f1.txt')
    os.remove(errors_file_path)

    main('f1', '1', '1', '10')

    errors = """\
1 0.596153846154
2 0.136752136752
3 0.00923076923077
4 4.06325614181e-17
5 1.18745157937e-16
6 2.01674114107e-16
7 1.12002254426e-16
8 4.92144386141e-17
9 2.63586228165e-17
10 7.8813157923e-19
"""

    with open(errors_file_path) as error_file:
        assert error_file.read() == errors
