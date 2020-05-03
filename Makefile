release:
ifndef VERSION
	$(error "Usage: make release VERSION=0.1.9")
endif
	git checkout master
	git pull
	echo "__version__ = '$(VERSION)'" > usum/__version__.py
	git add usum/__version__.py
	git commit -m "Set version to $(VERSION)"
	git push
	make dist
	twine upload dist/usum-$(VERSION)*
	git checkout develop
	git pull
	git rebase origin/master
	@echo "Create a new release version on: https://github.com/prihoda/usum/releases"

dist:
	python setup.py sdist bdist_wheel
