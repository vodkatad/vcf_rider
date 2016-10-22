#[derive(Debug)]
struct A {
    a: Vec<u8>,
    b: usize
}

fn main() {
    let names = vec![
        A { a: Vec::new(), b: 0 },
        A { a: Vec::new(), b: 1 },
        A { a: Vec::new(), b: 2 }
    ];

    let stuff = names
        .iter()
        .map(|x| x.b);
        //.fold(0, |acc, len| acc + len );

    //assert_eq!(total_bytes, 3);
    //use_names_for_something_else(names);
    for x in names {
        println!("{:?}", x);
    }
}

fn use_names_for_something_else(names: Vec<A>) {
    println!("{:?}", names);
}