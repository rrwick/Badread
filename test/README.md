# Badread tests

Badread comes with a few automated tests to help with development and spotting bugs.

To run the tests, execute this command from Badread's root directory:
```
python3 -m unittest
```

Or if you have [Coverage.py](https://coverage.readthedocs.io/en/coverage-4.5.1a/) installed, you can run the tests through it to get some code coverage stats:
```
coverage run -m unittest && coverage report -m
```

The tests aren't particularly fast, because of the stochastic nature of read simulation. E.g. some tests run hundreds of trials to make sure that simulated read identities stay within expected bounds. So it may take a couple minutes for the tests to complete. Sorry, I know that's not ideal for unit tests!
