## utils.R tests here

library(Flow)

context('check functions in utils.R')



test_that("file.name", {
    paths1 = c('path/path/long_path/usr/foo123.png', 'name.name/usr/bin/file_1.txt', 'usr/bin/file.named.like.this.txt')
    expect_match(file.name(paths1)[1], 'foo123.png')
    expect_match(file.name(paths1)[2], 'file_1.txt')
    expect_match(file.name(paths1)[3], 'file.named.like.this.txt')

})


test_that("file.dir", {
    paths1 = c('path/path/long_path/usr/foo123.png', 'name.name/usr/bin/file_1.txt', 'usr/bin/file.named.like.this.txt')
    expect_match(file.dir(paths1)[1], 'path/path/long_path/usr/')
    expect_match(file.dir(paths1)[2], 'name.name/usr/bin')
    expect_match(file.dir(paths1)[3], 'usr/bin/')

})


