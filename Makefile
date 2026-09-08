build:
	cd third_party/gsetl; cargo build --release

install:
	cp third_party/gsetl/target/release/gsetl /usr/bin/

bai:
	cd third_party/gsetl && cargo build --release
	cp third_party/gsetl/target/release/gsetl /usr/bin/

clean:
	rm -rf third_party/gsetl/target
