# Badread tests

Badread comes with a few automated tests to help with development and spotting bugs.

To run the tests, execute this command from Porechop's root directory:
```
python3 -m unittest
```

Or if you have [coverage]() installed, you can run the tests through it to get some code coverage stats:
```
coverage run --source badread -m unittest && coverage report -m
```

The tests aren't as particularly fast, because of the stochastic nature of read simulation. E.g. some tests run hundreds of trials to make sure that simulated read identities stay within expected bounds. So it may take about 5–10 minutes for the tests to complete. Sorry, I know that's not ideal for unit tests.
