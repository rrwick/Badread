# Contributing to Badread

Thanks for your interest in contributing to Badread! Please have a quick read of this page before making a change, so you know what to expect.



## Bug fixes

If you have discovered a bug, please report it via the GitHub [issue tracker](https://github.com/rrwick/Badread/issues).

Bug-fix [pull requests](https://github.com/rrwick/Badread/pulls) are always welcome! Please ensure that Badread's [automated tests](../test) still pass after you have made your changes. And whenever possible, you are encouraged to add new automated tests which cover your changes. Ideally, create a test which fails before your fix and passes after.



## New features

Badread's intended scope is to be a tool which generates synthetic read sets that can be used for testing other tools which take long reads as input. If you have added a new feature to Badread which stays within this scope, I would be interested in merging it in via a [pull request](https://github.com/rrwick/Badread/pulls)! Along with your new feature, please make any necessary updates to the [README.md](https://github.com/rrwick/Badread/blob/master/README.md) file to describe the feature and its usage. If possible, you are also encouraged to write [automated tests](../test) which cover your new feature.

In order to avoid [feature creep](https://en.wikipedia.org/wiki/Feature_creep), I may not accept pull requests which depart from Badread's scope. E.g. a feature for automatically designing new barcode sequences – certainly useful for some but outside the scope. However, you are always welcome to fork Badread and make such changes in your own forked copy.



## Other contribution notes

* Pull requests which improve performance (but do not change Badread's behaviour) are also welcome!
* Badread is not using feature branches. I will therefore merge pull requests into the master branch and make a new [release](https://github.com/rrwick/Badread/releases) when I'm confident that everything works as intended.
* Don't worry about changing Badread's version number in your pull request. I will bump the version when I make a new release.
* Adding new dependencies to Badread is discouraged. Please refrain from doing so unless it is justified.



## Code of conduct

### Our pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education, socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.

### Our standards

Examples of behaviour that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behaviour by participants include:

* The use of sexualised language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

### Our responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behaviour and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behaviour.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviours that they deem inappropriate, threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behaviour may be reported by contacting the project team at rrwick@gmail.com. All complaints will be reviewed and investigated and will result in a response that is deemed necessary and appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the Contributor Covenant version 1.4, available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

For answers to common questions about this code of conduct, see https://www.contributor-covenant.org/faq
