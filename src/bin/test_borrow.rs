fn main() {
    let mut a = vec!(1,2,3);
    {
        let b = &mut a;
        let c = b.pop().unwrap(); // so this changes the len of var a on the stack?
        println!("{}", c);
    }
    println!("{}", a.len());
}