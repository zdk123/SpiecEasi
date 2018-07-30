context("Graphs")

types <-  c("cluster", "erdos_renyi",  "hub",
            "scale_free", "block", "band")
D <- 10 ; e <- 10
gr_list <- lapply(types, make_graph, D=D, e=e, enforce=TRUE)

test_that("synthetic graphs contain expected graph attribute", {
    lapply(1:length(gr_list), function(i) {
            expect_match(attr(gr_list[[i]], "graph"), types[i])
            })

})

test_that("unrecognized graph method returns error", {
    expect_error(make_graph('foo', 10, 10), "*graph method foo not supported*")
})



tol <- 9
conds <- round(exp(seq(log(10), log(10000), length.out=10)))
prec_list <- lapply(gr_list, function(gr) {
                lapply(conds, function(kappa) {
                    graph2prec(gr, epsBin=tol, targetCondition=kappa)
                 })
             })
test_that("condition number of precision matrix is within eps of target condition", {
    lapply(1:length(prec_list), function(i) {
        lapply(1:length(prec_list[[i]]), function(j) {
            expect_equal(kappa(prec_list[[i]][[j]]), conds[j], tolerance=tol, scale=1)
        })
    })
})


test_that("graph is s3 class graph", {
    lapply(gr_list, function(gr) {
        expect_match(class(gr), 'graph')
    })
})


test_that("graph2prec returns error if input is not of class graph", {
    fakePrec <- matrix(0, 10,10)
    class(fakePrec) <- 'foo'
    expect_error(graph2prec(fakePrec), 'input is not a graph')
})

test_that("expected number of edges is properly enforced", {
    for (gr in gr_list) {
      expect_equal(edge_count(gr), e)
    }
})
