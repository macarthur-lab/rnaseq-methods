TAG1 = weisburd/olego:v1.1.7
TAG2 = weisburd/olego:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
