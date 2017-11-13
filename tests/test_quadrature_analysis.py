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
4 4.07201315935e-17
5 1.18745157937e-16
6 2.01674114107e-16
7 1.12089824602e-16
8 4.92144386141e-17
9 2.6271052641e-17
10 7.8813157923e-19
"""

    with open(errors_file_path) as error_file:
        assert error_file.read() == errors

def test_errors_produced_by_f6():
    errors_file_path = os.path.join('errors', 'f6.txt')
    os.remove(errors_file_path)

    main('f6', '1', '1', '20')

    errors = """\
1 1.0
2 0.359761434409
3 0.0542663665355
4 0.00362156280181
5 0.000225249571234
6 2.47356012228e-06
7 1.3475281515e-08
8 3.60270774865e-10
9 1.89443243235e-12
10 3.44732056933e-16
11 8.34897950384e-17
12 7.40635278566e-18
13 1.14461815778e-16
14 1.85158819642e-18
15 6.11024104817e-17
16 1.35502590738e-17
17 1.34913449039e-16
18 5.681009239e-17
19 1.09496192888e-16
20 5.46218517943e-17
"""

    with open(errors_file_path) as error_file:
        assert error_file.read() == errors
