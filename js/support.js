function debugThis(t, tk) {
    if (t >= tk) {
        debugger
    }
}
function printLog(t, el) {
    console.log(t);
    el.forEach(item => {
        console.log(` ${item}`);
    });
}
export {
    debugThis,
    printLog
}