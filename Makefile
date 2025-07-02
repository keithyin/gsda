build:
	cd gsetl; cargo build --release

install:
	cp gsetl/target/release/gsetl /usr/bin/

bai:
	cd gsetl && cargo build --release
	cp gsetl/target/release/gsetl /usr/bin/

clean:
	rm -rf gsetl/target