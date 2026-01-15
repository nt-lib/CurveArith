.PHONY: clean copy_remote test_remote 
REMOTE_DIR := /tmp/CurveArith

clean:
	find . -type f -name "*.sig" -delete

copy_remote:
	rsync -avz --delete ./ $(ssh):$(REMOTE_DIR)

test:
	cd tests && magma -n src/test_all.m

# Usage:
# make test_remote ssh="user@hostname"
test_remote: copy_remote
	ssh $(ssh) "cd $(REMOTE_DIR) && make test"



